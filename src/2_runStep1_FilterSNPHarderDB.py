from pathlib import Path
from bisect import bisect_right
import sqlite3
import argparse
import gzip
import re
import multiprocessing as mp

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output-prefix", required=True)
    parser.add_argument("-c", "--config", required=True)

    args = parser.parse_args()

    filter_vcf(
        input_vcf=args.input,
        output_prefix=args.output_prefix,
        config_path=args.config
    )