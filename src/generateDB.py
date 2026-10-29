import sqlite3
import argparse
import re

def create_snp_db(input_file, db_file, table_name, batch_size=100000):
    # Validate table name
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", table_name):
        raise ValueError(f"Invalid table name: {table_name}")

    conn = sqlite3.connect(db_file)
    cur = conn.cursor()

    cur.execute(f"""
        CREATE TABLE IF NOT EXISTS "{table_name}" (
            CHROM TEXT NOT NULL,
            START INTEGER NOT NULL,
            END INTEGER NOT NULL
        )
    """)

    # Replace existing records
    cur.execute(f'DELETE FROM "{table_name}"')

    batch = []

    with open(input_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            chrom, start, end = line.split()[:3]
            batch.append((chrom, int(start), int(end)))

            if len(batch) >= batch_size:
                cur.executemany(
                    f'INSERT INTO "{table_name}" VALUES (?, ?, ?)',
                    batch
                )
                batch.clear()

        if batch:
            cur.executemany(
                f'INSERT INTO "{table_name}" VALUES (?, ?, ?)',
                batch
            )

    conn.commit()

    index_name = f"idx_{table_name}_region"

    cur.execute(f"""
        CREATE INDEX IF NOT EXISTS "{index_name}"
        ON "{table_name}" (CHROM, START, END)
    """)

    conn.commit()
    conn.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-t", "--table", required=True)

    args = parser.parse_args()

    create_snp_db(
        input_file=args.input,
        db_file=args.output,
        table_name=args.table
    )