import re
import argparse
import logging
from pathlib import Path
from collections import defaultdict

import anndata as ad
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.io import mmwrite


META_COLUMNS = {
    "chr", "pos", "refbase", "altbase", "mode",
    "altCellType", "totalCellNumber", "totalReadsNumber",
    "isMultiAlleles", "isMultiCellType", "isCluster",
    "SNP", "INDEL", "otherFilter"
}


def parse_altbase(value):
    """A[24],C[1] -> ['A', 'C']"""
    return re.findall(r"([^,\[\]]+)\[\d+\]", str(value))


def parse_cell_entry(entry, allowed_alts):
    """
    CTCGAAAGTCAAAGCG-C-1
    -> (CTCGAAAGTCAAAGCG, C, 1)
    """
    entry = entry.strip()
    if not entry:
        return None

    for alt in sorted(allowed_alts, key=len, reverse=True):
        match = re.fullmatch(
            rf"(.+)-({re.escape(alt)})-(\d+)",
            entry
        )
        if match:
            return (
                match.group(1),
                match.group(2),
                int(match.group(3))
            )

    return None


def convert(input_file, output_dir, prefix="cellxsnp", binary=False):
    output_dir = Path(output_dir) / "result" / "mtx"
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Reading %s", input_file)
    df = pd.read_csv(input_file, sep="\t", dtype=str)

    celltype_columns = [
        c for c in df.columns if c not in META_COLUMNS
    ]

    observations = defaultdict(int)
    cell_types = defaultdict(set)
    variant_metadata = {}

    for _, record in df.iterrows():
        chrom = str(record["chr"])
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom

        pos = str(record["pos"])
        ref = str(record["refbase"])
        alts = parse_altbase(record["altbase"])

        if not alts:
            logging.warning(
                "Skipping variant %s:%s: invalid ALT",
                chrom, pos
            )
            continue

        for alt in alts:
            variant = f"{chrom}:{pos}:{ref}>{alt}"

            variant_metadata[variant] = {
                "chr": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "mode": record.get("mode"),
                "SNP": record.get("SNP"),
                "INDEL": record.get("INDEL"),
                "otherFilter": record.get("otherFilter")
            }

        for celltype in celltype_columns:
            value = record.get(celltype)

            if pd.isna(value) or ";" not in str(value):
                continue

            _, entries = str(value).split(";", 1)

            for entry in entries.split(":"):
                parsed = parse_cell_entry(entry, alts)

                if parsed is None:
                    if entry.strip():
                        logging.warning(
                            "Unparsed entry at %s:%s: %s",
                            chrom, pos, entry
                        )
                    continue

                barcode, alt, read_count = parsed
                variant = f"{chrom}:{pos}:{ref}>{alt}"

                observations[(barcode, variant)] += read_count
                cell_types[barcode].add(celltype)

    # Build cell and variant indices
    cells = sorted({
        cell for cell, _ in observations
    })
    variants = list(variant_metadata)

    cell_index = {c: i for i, c in enumerate(cells)}
    variant_index = {v: j for j, v in enumerate(variants)}

    rows, cols, values = [], [], []

    for (cell, variant), count in observations.items():
        rows.append(cell_index[cell])
        cols.append(variant_index[variant])
        values.append(1 if binary else count)

    matrix = coo_matrix(
        (values, (rows, cols)),
        shape=(len(cells), len(variants)),
        dtype="int32"
    ).tocsr()

    matrix.eliminate_zeros()

    # Cell metadata
    cells_df = pd.DataFrame({
        "cell_id": cells,
        "celltype": [
            ",".join(sorted(cell_types[c]))
            for c in cells
        ]
    })

    # Variant metadata
    variants_df = pd.DataFrame([
        {"variant_id": v, **variant_metadata[v]}
        for v in variants
    ])

    # 1. Save Matrix Market
    mtx_path = output_dir / f"{prefix}.mtx"
    mmwrite(str(mtx_path), matrix)

    # 2. Save cell and variant metadata
    cells_df.to_csv(
        output_dir / f"{prefix}.cells.tsv",
        sep="\t",
        index=False
    )

    variants_df.to_csv(
        output_dir / f"{prefix}.variants.tsv",
        sep="\t",
        index=False
    )

    # 3. Create AnnData
    adata = ad.AnnData(
        X=matrix,
        obs=cells_df.set_index("cell_id"),
        var=variants_df.set_index("variant_id")
    )

    adata.uns["matrix_type"] = (
        "binary" if binary else "alt_read_count"
    )

    # 4. Save H5AD
    h5ad_path = output_dir / f"{prefix}.h5ad"
    adata.write_h5ad(
        h5ad_path,
        compression="gzip"
    )

    # 5. Save Cell x SNP CSV
    # Sparse DataFrame avoids densifying the matrix in memory.
    csv_path = output_dir / f"{prefix}.csv"

    cellxsnp_df = pd.DataFrame.sparse.from_spmatrix(
        matrix,
        index=cells,
        columns=variants
    )

    cellxsnp_df.to_csv(
        csv_path,
        index_label="cell_id"
    )

    logging.info("Cells: %s", f"{len(cells):,}")
    logging.info("Variants: %s", f"{len(variants):,}")
    logging.info("Nonzero entries: %s", f"{matrix.nnz:,}")
    logging.info("Matrix shape: %s", matrix.shape)
    logging.info("Output directory: %s", output_dir)
    logging.info("MTX: %s", mtx_path)
    logging.info("H5AD: %s", h5ad_path)
    logging.info("CSV: %s", csv_path)

    return adata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert mutation table to Cell x SNP matrix"
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input mutation TSV file"
    )

    parser.add_argument(
        "--output",
        default="output",
        help="Base output directory (default: output)"
    )

    parser.add_argument(
        "--prefix",
        default="cellxsnp",
        help="Output file prefix"
    )

    parser.add_argument(
        "--binary",
        action="store_true",
        help="Store 0/1 mutation presence instead of ALT read counts"
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s"
    )

    convert(
        input_file=args.input,
        output_dir=args.output,
        prefix=args.prefix,
        binary=args.binary
    )
