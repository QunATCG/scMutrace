# import
import os
import sys
import pysam
import timeit
import argparse
from py_utils import make_dir_check
import multiprocessing as mp

'''
Option 1:
https://github.com/10XGenomics/subset-bam

Option 2:
There is also another method to extract bam with desired barcodes.
# filter.txt
CB:Z:GCCAAATTCACATACG-1
CB:Z:GGACATTGTGATGATA-1
CB:Z:TCAGGTACATTAGGCT-1

export BAM_FILE='possorted_genome_bam.bam'

# Save the header lines
samtools view -H $BAM_FILE > SAM_header

# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view $BAM_FILE | LC_ALL=C grep -F -f filter.txt > filtered_SAM_body

# Combine header and body
cat SAM_header filtered_SAM_body > filtered.sam

# Convert filtered.sam to BAM format
samtools view -b filtered.sam > filtered.bam

Option 3:
this script pysam
'''

# get argparse
def initialize_parser():
    parser = argparse.ArgumentParser(description='Script to exclude noisy reads based on barcode')
    parser.add_argument('--bam', type=str, help='BAM file to be analysed', required = True)
    parser.add_argument('--ref', type=str, help='Path to reference genome version. *fai must be available in the same directory as the reference genome file', required = True)
    parser.add_argument('--barcode', type=str, help='Barcode path', required = True)
    parser.add_argument('--barcodetag', type=str, help='Barcode tag', required = True)
    parser.add_argument('--contigFile', type=str, help='Contig information in tab-separated format: Chr\tstart\tend', required = True)
    parser.add_argument('--splitCheck', action='store_true', help='A boolean flag')
    parser.add_argument('--id', help='Prefix used to name output file.', required = True)
    parser.add_argument('--ncut', default = 10, help='Number of cut reads bases [Default: 10]',required=False,type = int)
    parser.add_argument('--nprocs', default = 1, help='Number of processes [Default: 1]',required=False,type = int)
    parser.add_argument('--outdir', default = '.', help='Out folder', required = False)
    return (parser)

def err_call_back(err):
    print(f'apply_async error:{str(err)}')

# get config data from fasta
def generate_config_range(refasta, contigSelfInPath):
    #basicChr = list(str(i) for i in range(1, 40)) + ["chr" + str(i) for i in range(1, 40)] + ['X', 'Y', 'chrX', 'chrY', 'M', 'chrM', 'm', 'chrm', 'MT', 'chrMT', 'mt', 'chrmt','x', 'y', 'chrx', 'chry']
    chrListSelfDefined = []
    with open(contigSelfInPath, "r") as f:
        for line in f:
            chrListSelfDefined.append(line.split("\t")[0].strip())
    refChrTotal = refasta.references
    configDic = {}
    chrList = []
    for refChr in refChrTotal:
        if refChr in chrListSelfDefined:
            sqlConfig = f"{refChr}_1_{refasta.get_reference_length(refChr)}"
            chrList.append(refChr)
            configDic[refChr] = sqlConfig
    return configDic, chrList

def mark_visium_duplicates(input_bam, output_bam, cb_tag="CB", ub_tag="UB", strip_cb_suffix=False, choose_best=True, tmp_dir=None, batch_size=100_000,):
    """
    Mark duplicate reads in Visium / Visium HD BAM using:
        CB + UB + genomic fragment position
    Paired-end:
        key = CB + UB + chr + fragment_start + fragment_end
    Single-end / abnormal pair:
        key = CB + UB + chr + start + end + strand
    Duplicate reads are marked with BAM flag 0x400.
    Reads are NOT removed.
    Memory-optimized version:
        - best_read is stored in SQLite instead of Python dict
        - duplicate query names are stored in SQLite instead of Python set
        - suitable for very large coordinate-sorted BAMs
    Parameters
    ----------
    input_bam : str
        Input coordinate-sorted BAM.
    output_bam : str
        Output BAM.
    cb_tag : str
        Cell/bin barcode tag, usually "CB".
    ub_tag : str
        Corrected UMI tag, usually "UB".
    strip_cb_suffix : bool
        Remove trailing "-1" etc. from CB if True.
    choose_best : bool
        If True, retain the fragment with:
            1. highest MAPQ
            2. highest total base quality
        If False, retain the first fragment encountered.
    tmp_dir : str or None
        Directory for temporary SQLite database.
        For HPC, ideally use fast local scratch storage.
    batch_size : int
        Number of records between SQLite commits.
    """
    # ============================================================
    # Temporary SQLite database
    # ============================================================
    if tmp_dir is None:
        tmp_dir = tempfile.gettempdir()
    os.makedirs(tmp_dir, exist_ok=True)
    fd, db_path = tempfile.mkstemp(prefix="visium_dedup_", suffix=".sqlite", dir=tmp_dir)
    os.close(fd)
    print(f"Temporary database: {db_path}")
    conn = sqlite3.connect(db_path)
    # Faster SQLite settings for temporary data
    conn.execute("PRAGMA journal_mode = OFF")
    conn.execute("PRAGMA synchronous = OFF")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute("PRAGMA locking_mode = EXCLUSIVE")

    # ============================================================
    # Table 1: best representative per molecule family
    # ============================================================
    print("\n[Pass 1/3] Finding representative reads...")
    conn.execute(
        """
        CREATE TABLE best_read (
            cb TEXT NOT NULL,
            ub TEXT NOT NULL,
            tid INTEGER NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            extra INTEGER NOT NULL,

            query_name TEXT NOT NULL,

            mapq INTEGER NOT NULL,
            baseq INTEGER NOT NULL,

            PRIMARY KEY (
                cb,
                ub,
                tid,
                start,
                end,
                extra
            )
        ) WITHOUT ROWID
        """
    )

    bam = pysam.AlignmentFile(input_bam, "rb")
    processed = 0
    usable = 0

    cursor = conn.cursor()

    for read in bam.fetch(until_eof=True):
        processed += 1
        if not valid_for_family(read):
            continue
        cb, ub = get_cb_ub(read, cb_tag, ub_tag, strip_cb_suffix)
        if cb is None:
            continue
        key = get_key(read, cb, ub)
        # paired read2
        if key is None:
            continue
        usable += 1
        mapq, baseq = get_score(read)
        (
            cb,
            ub,
            tid,
            start,
            end,
            extra,
        ) = key

        if choose_best:
            # Insert if absent.
            #
            # If family already exists:
            #   replace only when new MAPQ is better
            #   or MAPQ equal and baseQ is better.
            cursor.execute(
                """
                INSERT INTO best_read (
                    cb,
                    ub,
                    tid,
                    start,
                    end,
                    extra,
                    query_name,
                    mapq,
                    baseq
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)

                ON CONFLICT (
                    cb,
                    ub,
                    tid,
                    start,
                    end,
                    extra
                )

                DO UPDATE SET
                    query_name = excluded.query_name,
                    mapq = excluded.mapq,
                    baseq = excluded.baseq

                WHERE
                    excluded.mapq > best_read.mapq

                    OR (
                        excluded.mapq = best_read.mapq
                        AND
                        excluded.baseq > best_read.baseq
                    )
                """,
                (
                    cb,
                    ub,
                    tid,
                    start,
                    end,
                    extra,
                    read.query_name,
                    mapq,
                    baseq,
                ),
            )

        else:

            # Keep first encountered read
            cursor.execute(
                """
                INSERT OR IGNORE INTO best_read (
                    cb,
                    ub,
                    tid,
                    start,
                    end,
                    extra,
                    query_name,
                    mapq,
                    baseq
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    cb,
                    ub,
                    tid,
                    start,
                    end,
                    extra,
                    read.query_name,
                    mapq,
                    baseq,
                ),
            )
        if usable % batch_size == 0:
            conn.commit()
            print(
                f"\r  usable fragments: {usable:,}",
                end="",
                flush=True,
            )
    bam.close()
    conn.commit()
    print(
        f"\n  reads scanned:     {processed:,}"
        f"\n  fragments checked: {usable:,}"
    )

    n_families = conn.execute(
        "SELECT COUNT(*) FROM best_read"
    ).fetchone()[0]

    print(f"  molecule families: {n_families:,}")

    # ============================================================
    # Table 2: duplicate QNAMEs
    # ============================================================
    print("\n[Pass 2/3] Identifying duplicate query names...")
    conn.execute(
        """
        CREATE TABLE duplicate_names (
            query_name TEXT PRIMARY KEY
        ) WITHOUT ROWID
        """
    )
    bam = pysam.AlignmentFile(input_bam, "rb")
    checked = 0
    duplicate_fragments = 0
    select_best = conn.cursor()
    insert_dup = conn.cursor()
    for read in bam.fetch(until_eof=True):
        if not valid_for_family(read):
            continue
        cb, ub = get_cb_ub(read, cb_tag, ub_tag, strip_cb_suffix)
        if cb is None:
            continue
        key = get_key(read, cb, ub)
        if key is None:
            continue
        checked += 1

        (
            cb,
            ub,
            tid,
            start,
            end,
            extra,
        ) = key

        result = select_best.execute(
            """
            SELECT query_name
            FROM best_read
            WHERE
                cb = ?
                AND ub = ?
                AND tid = ?
                AND start = ?
                AND end = ?
                AND extra = ?
            """,
            (
                cb,
                ub,
                tid,
                start,
                end,
                extra,
            ),
        ).fetchone()

        if result is None:
            continue

        representative_name = result[0]

        if read.query_name != representative_name:
            insert_dup.execute(
                """
                INSERT OR IGNORE INTO duplicate_names
                VALUES (?)
                """,
                (read.query_name,),
            )
            duplicate_fragments += 1
        if checked % batch_size == 0:
            conn.commit()
            print(
                f"\r  fragments checked: {checked:,}",
                end="",
                flush=True,
            )
    bam.close()
    conn.commit()

    print(f"\n  duplicate fragments: {duplicate_fragments:,}")

    n_duplicate_names = conn.execute(
        "SELECT COUNT(*) FROM duplicate_names"
    ).fetchone()[0]

    print(f"  duplicate QNAMEs:     {n_duplicate_names:,}")

    # best_read is no longer needed.
    # Drop it to reduce DB size / lookup overhead.
    conn.execute("DROP TABLE best_read")
    conn.commit()

    # ============================================================
    # Pass 3:
    # mark BOTH mates using query_name
    # ============================================================

    print("\n[Pass 3/3] Writing duplicate-marked BAM...")

    bam = pysam.AlignmentFile(input_bam, "rb")
    out = pysam.AlignmentFile(
        output_bam,
        "wb",
        template=bam,
    )
    lookup_dup = conn.cursor()
    total_reads = 0
    duplicate_reads = 0
    for read in bam.fetch(until_eof=True):
        total_reads += 1
        result = lookup_dup.execute(
            """
            SELECT 1
            FROM duplicate_names
            WHERE query_name = ?
            """,
            (read.query_name,),
        ).fetchone()

        if result is not None:
            read.is_duplicate = True
            duplicate_reads += 1
        out.write(read)
        if total_reads % 1_000_000 == 0:
            print(
                f"\r  reads written: {total_reads:,}",
                end="",
                flush=True,
            )
    bam.close()
    out.close()
    conn.close()

    # ============================================================
    # Index
    # ============================================================

    print("Indexing output BAM...")
    pysam.index(output_bam)

    # Remove temporary database
    try:
        os.remove(db_path)
    except OSError:
        pass

    # ============================================================
    # Summary
    # ============================================================

    print("\nDone.")
    print(f"Total reads:      {total_reads:,}")
    print(f"Duplicate reads:  {duplicate_reads:,}")

    if total_reads > 0:
        print(
            f"Duplicate rate:   "
            f"{duplicate_reads / total_reads:.2%}"
        )
    else:
        print("Duplicate rate:   NA")

    print(f"Output:           {output_bam}")

# split bam based on barcode
def split_bam_basedList(bam_Path, barcodeSubList, ChrList, configValue, subam_path, seqTag, nameOut, ncut, barcodetag):
    bam = pysam.AlignmentFile(bam_Path, 'rb')
    count_num = 0
    chrtmp = str(configValue.split("_")[0])
    startmp = str(configValue.split("_")[1])
    endtmp = str(configValue.split("_")[2])
    # create and open small pysam files
    outBamPath = os.path.join(subam_path, f"{nameOut}_{chrtmp}_{startmp}_{endtmp}.tmp.{seqTag}.bam")
    Outfile = pysam.AlignmentFile(outBamPath, "wb", template=bam)
    print("###########################################")
    print("Split " + seqTag)
    print("Creating tmp bam files in " + outBamPath)
    print("###########################################")
    for read in bam.fetch(contig=chrtmp, start=int(startmp), stop=int(endtmp)):
        count_num += 1
        if count_num % 500000 == 0:
            print("Processed " + str(count_num) + " reads")
        if read.has_tag(barcodetag):
            cb_read = read.get_tag(barcodetag).split('-')[0]
            chr_read = read.reference_name
            if cb_read in barcodeSubList and chr_read in ChrList:
                # cut each read's first left and right bases
                CIGAR_read = read.cigartuples
                # check if CIGAR exists
                if CIGAR_read and len(CIGAR_read) > 1:
                    cut_l = CIGAR_read[0][1] + ncut if CIGAR_read[0][0] == 4 else ncut
                    cut_r = CIGAR_read[-1][1] + ncut if CIGAR_read[-1][0] == 4 else ncut
                else:
                    cut_l = ncut
                    cut_r = ncut
                base_quality = read.query_qualities
                index_cut_bq = list(range(0, cut_l)) + list(range(-1, -cut_r - 1, -1))
                if base_quality is not None:
                    for index_q in index_cut_bq:
                        base_quality[index_q] = 0
                    read.query_qualities = base_quality
                Outfile.write(read)
    Outfile.close()
    print("Indexing bam --- " + outBamPath)
    pysam.index(outBamPath)
    print("Split bam ###---### Done!")
    bam.close()

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    bam_Path = args.bam
    ref_Path = args.ref
    barcodePath = args.barcode
    barcodetag = args.barcodetag
    contigSelfInPath = args.contigFile
    split_check = args.splitCheck
    nameSplitBamTag = args.id
    seqTag = "SMT"
    ncut = args.ncut
    COREUSE = args.nprocs
    outDirPath = args.outdir
    # make dir
    make_dir_check(outDirPath)
    # extract barcode
    barcodeListCell = []
    barcodeListFile = open(outDirPath + "/barcodeList.txt", "w+")
    totalBarcode = open(barcodePath, "r")
    for line in totalBarcode:
        barcodeTmp = line.split()[0].strip().split("-")[0]
        if barcodeTmp not in barcodeListCell:
            barcodeListFile.writelines(barcodeTmp + "\n")
            barcodeListCell.append(barcodeTmp)
    barcodeListFile.close()
    # generate_config_range
    configDic, chrList = generate_config_range(pysam.FastaFile(ref_Path), contigSelfInPath)
    if len(chrList) == 0:
        print("No chromosome in the reference genome")
        sys.exit("No chr in the contig file")
    # extract bam by cellbarcode
    print(split_check)
    if split_check:
        if COREUSE > 1:
            print('Parent process %s.' % os.getpid())
            p = mp.get_context("fork").Pool(COREUSE)
            print("Your cores is " + str(COREUSE))
            results = []
            for value in configDic.values():
                p.apply_async(split_bam_basedList, args = (bam_Path, barcodeListCell, chrList, value, outDirPath, seqTag, nameSplitBamTag, ncut, barcodetag,), error_callback = err_call_back)
            print('Waiting for all subprocesses done...')
            p.close()
            p.join()
            for result in results:
                print(result.get())
            print('All subprocesses done.')
        else:
            for value in configDic.values():
                split_bam_basedList(bam_Path = bam_Path, barcodeSubList = barcodeListCell, ChrList = chrList, configValue = value, subam_path = outDirPath, seqTag = seqTag, nameOut = nameSplitBamTag, ncut=ncut, barcodetag=barcodetag)

if __name__=='__main__':
    start = timeit.default_timer()
    #19_1_59128983
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')