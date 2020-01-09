#!/usr/bin/env python3
import os
import re
import sys
import tempfile
import argparse
import logging
import subprocess


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("hisat2")


def run(cmd):
    logger.debug(f"Running command: {cmd}")
    subprocess.run(cmd, check=True, shell=True)

def hisat2(ubam, readgroup, bam, index, rna_strandness, summary_file, threads, temp_dir):
    with tempfile.TemporaryDirectory(dir=temp_dir) as tmpdir:
        logger.debug(f"temp directory = {tmpdir}")
        name = re.sub(".unmapped.bam$", "", os.path.basename(ubam))
        r1 = os.path.join(tmpdir, f"{name}_1.fq.gz")
        r2 = os.path.join(tmpdir, f"{name}_2.fq.gz")
        os.mkfifo(r1)
        os.mkfifo(r2)
        run(f"samtools fastq -1 {r1} -2 /dev/null -0 /dev/null -s /dev/null -O -n -F 0x900 {ubam} 2>/dev/null &")
        run(f"samtools fastq -1 /dev/null -2 {r2} -0 /dev/null -s /dev/null -O -n -F 0x900 {ubam} 2>/dev/null &")

        cmd = f"""hisat2 -x {index} -1 {r1} -2 {r2}"""
        if threads:
            cmd += f" -p {threads}"
        if rna_strandness:
            cmd += f" --rna-strandness {rna_strandness}"
        if summary_file:
            cmd += f" --summary-file {summary_file} --new-summary"

        if readgroup:
            cmd += f" --rg-id {readgroup['ID']}"
            keys = ["SM", "LB", "PL", "PU", "CN", "DT", "PM", "DS"]
            for key in keys:
                if key in readgroup:
                    value = readgroup[key]
                    if key == "DS":
                        # hisat2 hangs when there are pipes in the DS value.
                        # This is because it is not quoting it correctly when
                        # passing to subprograms, it also won't handle other
                        # special characters or spaces in file names etc.
                        # Because of the way hisat2 processing the command line
                        # it is difficult to quote correctly without quite a
                        # lot of changes to hisat2.
                        value = value.replace("|", "!")
                    cmd += f" --rg '{key}:{value}'"

        cmd += f" 2>/dev/null | samtools sort -o {bam}"
        run(cmd)
        cmd = f"samtools index -b {bam}"
        run(cmd)

def readgroup_str_to_dict(rg):
    return dict(x.split(":", 1) for x in rg.split("\t")[1:])

def bam_readgroups(fn):
    p = subprocess.run(
        f"samtools view -H {fn}",
        check=True,
        capture_output=True,
        shell=True)
    output = p.stdout.decode("ascii").splitlines()
    return [readgroup_str_to_dict(x) for x in output if x.startswith("@RG")]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", dest="bam", required=True)
    parser.add_argument("-i", "--index", required=True)
    parser.add_argument("-s", "--rna-strandness", choices=["F", "R", "FR", "RF"])
    parser.add_argument("-S", "--summary-file")
    parser.add_argument("-p", dest="threads")
    parser.add_argument("-T", "--temp-dir", default=tempfile.gettempdir())
    parser.add_argument("ubam")
    args = parser.parse_args()

    if not args.ubam.endswith(".unmapped.bam"):
        logger.error(f"uBAM file name must end with .unmapped.bam was {args.ubam}")
        sys.exit(1)
    logger.info(f"Input uBAM: {args.ubam}")

    # pysam doesn't seem to install on python 3.8 and is a pita to install on
    # alpine. Should it just call out to samtools to get the readgroup?

    readgroups = bam_readgroups(args.ubam)

    if len(readgroups) == 0:
        logger.warn(f"uBAM does not have a read group, no readgroup will be in output")

    if len(readgroups) > 1:
        logger.error(f"uBAM has multiple read groups")
        sys.exit(1)

    readgroup = readgroups[0] if len(readgroups) == 1 else {}

    hisat2(
        os.path.abspath(args.ubam),
        readgroup,
        os.path.abspath(args.bam),
        os.path.abspath(args.index),
        args.rna_strandness,
        args.summary_file,
        args.threads,
        args.temp_dir
    )


if __name__ == "__main__":
    main()

