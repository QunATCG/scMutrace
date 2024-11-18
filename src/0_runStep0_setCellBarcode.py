import pysam
import timeit
import argparse

def initialize_parser():
    parser = argparse.ArgumentParser(description='Script to set cell barcode')
    parser.add_argument('--bam', type=str, default='', help='BAM file to be analysed', required = True)
    parser.add_argument('--barcode', type=str, default='', help='Barcode file', required = True)
    parser.add_argument('--outpath', default = '.', help='Out folder', required = True)
    return (parser)

def add_cell_barcode(input_bam, output_bam, barcode):
    """
    Add a cell barcode to each read in a BAM file.

    Parameters:
    input_bam (str): Path to the input BAM file.
    output_bam (str): Path to the output BAM file.
    barcode (str): The cell barcode to add to each read.
    """
    # Open the input BAM file
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        # Create a new BAM file to write to
        with pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
            # Iterate over each read in the input BAM file
            for read in infile:
                # Set the cell barcode tag (CB) to the specified barcode
                read.set_tag("CB", barcode, value_type='Z')  # 'Z' is the type for a string tag

                # Write the modified read to the output BAM file
                outfile.write(read)

    print(f"Cell barcode '{barcode}' added to all reads in '{output_bam}'.")

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    inbam = args.bam
    barcode = args.barcode
    outbam = args.outpath
    add_cell_barcode(inbam, outbam, barcode)

if __name__ == "__main__":
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')