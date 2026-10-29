# import
import os
import sys
import pysam
import timeit
import argparse
import gc
import sqlite3
import tempfile
import traceback
import contextlib
from collections import defaultdict
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
    parser = argparse.ArgumentParser(description='Script to exclude noisy and duplicate reads based on barcode')
    parser.add_argument('--bam', type=str, help='BAM file to be analysed', required = True)
    parser.add_argument('--ref', type=str, help='Path to reference genome version. *fai must be available in the same directory as the reference genome file', required = True)
    parser.add_argument('--barcode', type=str, help='Barcode path', required = False)
    parser.add_argument('--barcodetag', type=str, help='Barcode tag', required = True)
    parser.add_argument('--umitag', type=str, help='UMI tag', required = True)
    parser.add_argument('--contigFile', type=str, help='Contig information in tab-separated format: Chr\tstart\tend', required = True)
    parser.add_argument('--splitCheck', action='store_true', help='A boolean flag')
    parser.add_argument('--id', help='Prefix used to name output file.', required = True)
    parser.add_argument('--ncut', default = 10, help='Number of cut reads bases [Default: 10]',required=False,type = int)
    parser.add_argument('--nprocs', default = 1, help='Number of processes [Default: 1]',required=False,type = int)
    parser.add_argument('--outdir', default = '.', help='Out folder', required = False)
    return (parser)

def finished_callback(result):
    print(
        f"[DONE] {result['region']} | "
        f"scanedReads={result['reads_scanned']:,} | "
        f"log={result['log']}",
        flush=True,
    )

def err_call_back(error):
    print(
        f"[ERROR] {error}",
        flush=True
    )


def get_cb_ub(read, cb_tag, ub_tag, strip_cb_suffix = False):
    if not read.has_tag(cb_tag):
        return None, None
    if not read.has_tag(ub_tag):
        return None, None
    cb = read.get_tag(cb_tag)
    ub = read.get_tag(ub_tag)
    if strip_cb_suffix:
        cb = cb.rsplit("-", 1)[0]
    return cb, ub

def valid_for_family(read):
    """
    Reads used to define duplicate families.
    """
    if read.is_unmapped:
        return False
    # Do not let secondary/supplementary alignments
    # define molecule families.
    if read.is_secondary:
        return False
    if read.is_supplementary:
        return False
    return True

def get_key(read, cb, ub):
    """
    Return a normalized SQLite key.
    Also returns False if this read should not be used
    to define a paired-end molecule (e.g. read2).
    """
    # --------------------------------------------------------
    # Normal paired-end fragment
    # --------------------------------------------------------
    if (read.is_paired and not read.mate_is_unmapped and read.reference_id == read.next_reference_id and read.template_length != 0):
        # Use read1 only to define fragment family
        if not read.is_read1:
            return None
        fragment_start = min(read.reference_start, read.next_reference_start)
        fragment_end = (fragment_start + abs(read.template_length))
        # final element pair_type = 0
        return (cb, ub, read.reference_id, fragment_start, fragment_end, 0)
    # --------------------------------------------------------
    # Single-end / abnormal pair
    # --------------------------------------------------------
    else:
        return (cb, ub, read.reference_id, read.reference_start, read.reference_end, 1 if read.is_reverse else 0)

def get_score(read):
    if read.query_qualities is None:
        base_quality_score = 0
    else:
        base_quality_score = sum(read.query_qualities)
    return (int(read.mapping_quality), int(base_quality_score))

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

# split bam based on barcode
def split_bam_basedList_v0(bam_Path, barcodeSubList, ChrList, configValue, subam_path, seqTag, nameOut, ncut, barcodetag):
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
def split_bam_basedList(bam_Path, barcodeSubList, ChrList, configValue, subam_path, seqTag, nameOut, ncut, barcodetag, ub_tag, tmp_dir, strip_cb_suffix=False, choose_best=True):
    """
    split bam based on barcode and chromosome list, and cut the first ncut bases of each read's left and right ends
    Parameters
    ----------
    bam_Path : str
        Input coordinate-sorted BAM.
    barcodeSubList : list|None
        List of barcodes to be included in the split BAM.
    ChrList : list
        List of chromosomes to be included in the split BAM.
    configValue : str
        Configuration value in the format "chr_start_end" to specify the region to be included in the split BAM.
    subam_path : str
        Path to the directory where the split BAM files will be saved.
    seqTag : str
        Sequence tag to be included in the output BAM file name.
    nameOut : str
        Prefix used to name output file.
    ncut : int
        Number of bases to be cut from the left and right ends of each read.
    barcodetag : str
        Barcode tag to be used for filtering reads.
    """
    chrtmp, startmp, endtmp = configValue.split("_", 2)
    startmp = int(startmp)
    endtmp = int(endtmp)
    # create and open small pysam files
    outBamPath = os.path.join(subam_path, f"{nameOut}_{chrtmp}_{startmp}_{endtmp}.tmp.{seqTag}.bam")
    rmdupBamPath = outBamPath.replace(".tmp.",".rmdup.")
    log_dir = os.path.join(subam_path, "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, f"step0_{nameOut}_{chrtmp}_{startmp}_{endtmp}.{seqTag}.splitrmdup.log")
    with open(log_path, "w", buffering=1) as log_file:
        with (contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file)):
            start_time = timeit.default_timer()
            try:
                print("=" * 80, flush=True)
                print("scMutrace mutation calling")
                print(f"PID:          {os.getpid()}", flush=True)
                print(f"Region:       {chrtmp}:{startmp}-{endtmp}", flush=True)
                print(f"Seq tag:      {seqTag}", flush=True)
                print(f"Input BAM:    {bam_Path}", flush=True)
                print(f"Temp BAM:     {outBamPath}", flush=True)
                print(f"Output BAM:   {rmdupBamPath}", flush=True)
                print("=" * 80, flush=True)
                # =================================================
                # Step 1: split BAM
                # =================================================
                print("\n[Step 1/2] Splitting BAM...", flush=True)
                bam = pysam.AlignmentFile(bam_Path, 'rb')
                Outfile = pysam.AlignmentFile(outBamPath, "wb", template=bam)
                count_num = 0
                if barcodeSubList is not None:
                    barcode_set = set(barcodeSubList)
                try:
                    for read in bam.fetch(contig=chrtmp, start=startmp, stop=endtmp):
                        count_num += 1
                        if count_num % 5000000 == 0:
                            print(f"Processed {count_num:,} reads", flush=True)
                        if not read.has_tag(barcodetag):
                            continue
                        cb_read = read.get_tag(barcodetag).split('-')[0]
                        if barcodeSubList is not None:
                            if cb_read not in barcode_set:
                                continue
                        if read.reference_name not in ChrList:
                            continue
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
                        if base_quality is not None:
                            read_length = len(base_quality)
                            left_end = min(cut_l, read_length)
                            right_start = max(0, read_length - cut_r)
                            # Set left terminal qualities to 0
                            for i in range(left_end):
                                base_quality[i] = 0
                            # Set right terminal qualities to 0
                            for i in range(right_start, read_length):
                                base_quality[i] = 0
                            read.query_qualities = (base_quality)
                        Outfile.write(read)
                finally:
                    Outfile.close()
                    bam.close()
                print(f"\nSplit BAM finished.", flush=True)
                print(f"  Reads scanned:  {count_num:,}",flush=True)
                # ------------------------------------------
                # Index temporary BAM
                # ------------------------------------------
                print("Indexing temporary BAM...", flush=True)
                pysam.index(outBamPath)
                # =================================================
                # Step 2: mark duplicates
                # =================================================
                print("\n[Step 2/2] Marking duplicates...", flush=True)
                mark_visium_duplicates(outBamPath, rmdupBamPath, barcodetag, ub_tag, strip_cb_suffix, choose_best, tmp_dir, batch_size=100_000,)
                # remove tmp bam file
                os.remove(outBamPath)
                os.remove(f"{outBamPath}.bai")
                elapsed = timeit.default_timer() - start_time
                print("\n" + "=" * 80)
                print("TASK FINISHED SUCCESSFULLY")
                print(
                    f"Region:       "
                    f"{chrtmp}:{startmp}-{endtmp}")
                print(
                    f"Elapsed:      "
                    f"{elapsed / 60:.2f} min")
                print(
                    f"Output BAM:   "
                    f"{rmdupBamPath}")
                print(
                    f"Log:          "
                    f"{log_path}")
                print("=" * 80, flush=True)
                # return value to apply_async callback
                return {
                    "region": configValue,
                    "log": log_path,
                    "reads_scanned": count_num,
                }
            except Exception:
                print("\n" + "=" * 80, flush=True)
                print("ERROR", flush=True)
                print(f"Region: {configValue}", flush=True)
                print("=" * 80, flush=True)
                traceback.print_exc(file=log_file)
                log_file.flush()
                # Important:
                # re-raise so apply_async(error_callback=...)
                # still knows this worker failed.
                raise

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    bam_Path = args.bam
    ref_Path = args.ref
    barcodePath = args.barcode
    barcodetag = args.barcodetag
    umitag = args.umitag
    contigSelfInPath = args.contigFile
    split_check = args.splitCheck
    nameSplitBamTag = args.id
    seqTag = "SMT"
    ncut = args.ncut
    COREUSE = args.nprocs
    outDirPath = args.outdir
    make_dir_check(outDirPath)
    tmp_dir = os.path.join(outDirPath, "tmp")
    make_dir_check(tmp_dir)
    # extract barcode
    barcodeListCell = set()
    if barcodePath is None:
        barcodeListCell = None
    else:
        with open(barcodePath, "r") as fin, open(f"{outDirPath}/barcodeList.txt", "w") as fout:
            for line in fin:
                barcodeTmp = line.split()[0].strip().split("-")[0]
                if barcodeTmp not in barcodeListCell:
                    fout.write(barcodeTmp + "\n")
                    barcodeListCell.add(barcodeTmp)
    configDic, chrList = generate_config_range(pysam.FastaFile(ref_Path), contigSelfInPath)
    if len(chrList) == 0:
        print("No chromosome in the reference genome")
        sys.exit("No chr in the contig file")
    if split_check:
        if COREUSE > 1:
            print('Parent process %s.' % os.getpid())
            num_tasks = len(configDic)
            cpu_count = mp.cpu_count()
            COREUSE = min(COREUSE, num_tasks, cpu_count)
            COREUSE = max(1, COREUSE)
            p = mp.get_context("fork").Pool(COREUSE)
            print("Your cores is " + str(COREUSE))
            results = []
            for value in configDic.values():
                result = p.apply_async(split_bam_basedList, args = (bam_Path, barcodeListCell, chrList, value, outDirPath, seqTag, nameSplitBamTag, ncut, barcodetag, umitag, tmp_dir, False, True,), callback=finished_callback, error_callback = err_call_back)
                results.append(result)
            p.close()
            print('Waiting for all subprocesses done...')
            p.join()
            for result in results:
                result.get()
            print('All subprocesses done.')
        else:
            for value in configDic.values():
                split_bam_basedList(bam_Path = bam_Path, barcodeSubList = barcodeListCell, ChrList = chrList, configValue = value, subam_path = outDirPath, seqTag = seqTag, nameOut = nameSplitBamTag, ncut=ncut, barcodetag=barcodetag, ub_tag=umitag, tmp_dir=tmp_dir, strip_cb_suffix=False, choose_best=True)
        print("Merging BAM files...")
        bam_list = []
        for value in configDic.values():
            chrtmp, startmp, endtmp = value.split("_", 2)
            bam_file = os.path.join(outDirPath, f"{nameSplitBamTag}_{chrtmp}_{startmp}_{endtmp}.rmdup.{seqTag}.bam")
            if not os.path.exists(bam_file):
                raise FileNotFoundError(f"Missing BAM file: {bam_file}")
            bam_list.append(bam_file)
        merged_bam = os.path.join(outDirPath, f"{nameSplitBamTag}.rmdup.{seqTag}.bam")
        pysam.merge("-@", str(COREUSE), "-f", "-c", merged_bam, *bam_list)
        pysam.index("-@", str(COREUSE), merged_bam)
        print("Merge finished.")
        print("Sorting BAM files...")
        sorted_bam = os.path.join(outDirPath, f"{nameSplitBamTag}.rmdup.sorted.{seqTag}.bam")
        pysam.sort("-@", str(COREUSE), "-o", sorted_bam, merged_bam)
        pysam.index("-@", str(COREUSE), sorted_bam)
        print("Sort finished.")
        os.remove(merged_bam)
        os.remove(f"{merged_bam}.bai")
        for bam in bam_list:
            os.remove(bam)
            os.remove(f"{bam}.bai")

if __name__=='__main__':
    start = timeit.default_timer()
    #19_1_59128983
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')