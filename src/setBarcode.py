import pysam
import timeit
import argparse
import os

def initialize_parser():
    parser = argparse.ArgumentParser(description='Add CB tag to BAM reads with buffered write/SMARTSEQ')
    parser.add_argument('--bam', type=str, help='Input BAM file', required=True)
    parser.add_argument('--outbam', type=str, help='Output BAM file', required=True)
    parser.add_argument('--buffer_size', type=int, default=500000, help='Number of reads to buffer before writing')
    return parser

def write_reads(outfile, buffer):
    for read in buffer:
        outfile.write(read)

def add_cell_barcode(input_bam, output_bam, buffer_size):
    try:
        with pysam.AlignmentFile(input_bam, "rb") as infile:
            with pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
                buffer = []
                for i, read in enumerate(infile, 1):
                    barcode = read.query_name.split(".")[0]
                    read.set_tag("CB", barcode, value_type='Z')
                    buffer.append(read)

                    if len(buffer) >= buffer_size:
                        write_reads(outfile, buffer)
                        buffer = []

                    if i % 1000000 == 0:
                        print(f"{i} reads processed...")

                if buffer:
                    write_reads(outfile, buffer)

        pysam.index(output_bam)
    except Exception as e:
        print(f"Error: {e}")

def main():
    parser = initialize_parser()
    args = parser.parse_args()

    inbam = args.bam
    outbam = args.outbam
    buffer_size = args.buffer_size

    add_cell_barcode(inbam, outbam, buffer_size)

if __name__ == "__main__":
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print(f'Computation time: {round(stop - start)} seconds')
