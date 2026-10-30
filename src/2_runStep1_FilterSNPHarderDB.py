from pathlib import Path
from bisect import bisect_right
import sqlite3
import argparse
import gzip
import re
import logging
import time

def normalize_chrom(chrom):
    chrom = str(chrom).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom.upper() in ("M", "MT"):
        return "chrM"
    if chrom.upper() in ("X", "Y"):
        return "chr" + chrom.upper()
    if chrom.isdigit():
        return "chr" + str(int(chrom))
    return "chr" + chrom

def open_vcf(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

class RegionDB_v0:
    # need big mem
    def __init__(self, name, db_path, table_name, action):
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", table_name):
            raise ValueError(f"Invalid table name: {table_name}")
        self.name = name
        self.action = action
        db_path = Path(db_path).resolve()
        self.conn = sqlite3.connect(
            f"file:{db_path}?mode=ro",
            uri=True
        )
        self.cur = self.conn.cursor()
        self.sql = f"""
            SELECT 1
            FROM "{table_name}"
            WHERE CHROM = ?
              AND START <= ?
              AND END >= ?
            LIMIT 1
        """
    def contains(self, chrom, pos):
        chrom = normalize_chrom(chrom)
        self.cur.execute(self.sql, (chrom, pos, pos))
        return self.cur.fetchone() is not None
    def close(self):
        self.conn.close()

class RegionDB:
    def __init__(self, name, db_path, table_name, action):
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", table_name):
            raise ValueError(f"Invalid table name: {table_name}")

        self.name = name
        self.action = action
        self.table_name = table_name

        db_path = Path(db_path).resolve()
        self.conn = sqlite3.connect(
            f"file:{db_path}?mode=ro", uri=True
        )

        self.current_chrom = None
        self.starts = []
        self.ends = []

    def load_chrom(self, chrom):
        chrom = normalize_chrom(chrom)

        sql = f"""
            SELECT START, END
            FROM "{self.table_name}"
            WHERE CHROM = ?
            ORDER BY START, END
        """

        starts, ends = [], []

        for start, end in self.conn.execute(sql, (chrom,)):
            start, end = int(start), int(end)

            if not starts or start > ends[-1] + 1:
                starts.append(start)
                ends.append(end)
            else:
                ends[-1] = max(ends[-1], end)

        self.current_chrom = chrom
        self.starts = starts
        self.ends = ends

    def contains(self, chrom, pos):
        chrom = normalize_chrom(chrom)

        if chrom != self.current_chrom:
            self.load_chrom(chrom)

        idx = bisect_right(self.starts, pos) - 1
        return idx >= 0 and pos <= self.ends[idx]

    def close(self):
        self.conn.close()


def load_databases(config_path):
    databases = []
    config_path = Path(config_path).resolve()
    with open(config_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            name, db_path, table_name, action = line.split()
            if action not in {"filter", "germline", "annotate"}:
                raise ValueError(f"Invalid action: {action}")
            db_path = Path(db_path)
            if not db_path.is_absolute():
                db_path = config_path.parent / db_path
            databases.append(
                RegionDB(name, db_path, table_name, action)
            )
    return databases

def annotate_info(hits):
    if not hits:
        return "None"
    annotation = ",".join(hits)
    return annotation

def filter_vcf(input_vcf, output_prefix, config_path):
    databases = load_databases(config_path)
    output_files = {
        "somatic": f"{output_prefix}.somatic_candidates.vcf",
        "germline": f"{output_prefix}.germline_candidates.vcf",
        "filtered": f"{output_prefix}.filtered.vcf"
    }
    counts = {
        "total": 0,
        "somatic": 0,
        "germline": 0,
        "filtered": 0
    }
    try:
        with open_vcf(input_vcf) as fin, \
             open(output_files["somatic"], "w") as fs, \
             open(output_files["germline"], "w") as fg, \
             open(output_files["filtered"], "w") as ff:
            outputs = [fs, fg, ff]
            for line in fin:
                if line.startswith("##"):
                    for out in outputs:
                        out.write(line)
                    continue
                if line.startswith("#CHROM"):
                    info_header = (
                        '##INFO=<ID=REGION_DB,Number=.,Type=String,'
                        'Description="Matched region databases">\n'
                    )
                    for out in outputs:
                        out.write(info_header)
                        out.write(line)
                    continue
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                chrom = fields[0]
                pos = int(fields[1])
                counts["total"] += 1
                hits = []
                actions = set()
                for db in databases:
                    if db.contains(chrom, pos):
                        hits.append(db.name)
                        actions.add(db.action)
                ann_info = annotate_info(hits)
                output_line = line.strip() + "\t" + ann_info + "\n"
                # Priority: filter > germline > somatic
                if "filter" in actions:
                    ff.write(output_line)
                    counts["filtered"] += 1

                elif "germline" in actions:
                    fg.write(output_line)
                    counts["germline"] += 1

                else:
                    fs.write(output_line)
                    counts["somatic"] += 1

    finally:
        for db in databases:
            db.close()

    print("\nFiltering summary")
    print("-----------------------------")
    for key, value in counts.items():
        print(f"{key:12s}: {value:,}")
    return counts

def setup_logging(log_dir):
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"step3_filterDB.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, encoding="utf-8"),
            logging.StreamHandler()
        ],
        force=True
    )

    return log_file


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input-dir", required=True)
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("-l", "--log-dir", default=None)

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    log_dir = (
        Path(args.log_dir)
        if args.log_dir
        else output_dir / "logs"
    )

    log_file = setup_logging(log_dir)

    vcf_files = sorted(input_dir.glob("*_filter.vcf"))

    logging.info("VCF filtering started")
    logging.info("Input directory: %s", input_dir)
    logging.info("Output directory: %s", output_dir)
    logging.info("Database config: %s", args.config)
    logging.info("Number of VCF files: %d", len(vcf_files))

    if not vcf_files:
        logging.error("No VCF files found in %s", input_dir)
        raise SystemExit(1)

    total_start = time.perf_counter()
    failed_files = 0

    for vcf_file in vcf_files:

        logging.info("=" * 60)
        logging.info("Processing: %s", vcf_file.name)

        start = time.perf_counter()

        output_prefix = output_dir / vcf_file.stem

        try:
            counts = filter_vcf(
                input_vcf=vcf_file,
                output_prefix=output_prefix,
                config_path=args.config
            )

            elapsed = time.perf_counter() - start


            logging.info(
                "SUCCESS: %s | Total=%d | Somatic=%d | "
                "Germline=%d | Filtered=%d | Time=%.2fs",
                vcf_file.name,
                counts["total"],
                counts["somatic"],
                counts["germline"],
                counts["filtered"],
                elapsed
            )

        except Exception:
            failed_files += 1
            elapsed = time.perf_counter() - start

            logging.exception(
                "FAILED: %s after %.2f seconds",
                vcf_file.name,
                elapsed
            )


    total_elapsed = time.perf_counter() - total_start

    logging.info("=" * 60)
    logging.info("All VCF files processed")
    logging.info("Total files: %d", len(vcf_files))
    logging.info("Successful: %d", len(vcf_files) - failed_files)
    logging.info("Failed: %d", failed_files)
    logging.info("Total elapsed: %.2f seconds", total_elapsed)
    logging.info("Log file: %s", log_file)

    if failed_files:
        raise SystemExit(1)