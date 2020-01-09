#!/usr/bin/env python3
import os
import re
import sys
import tempfile
import argparse
import logging
import subprocess


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("hisat2")


def run(cmd):
    logger.debug(f"Running command: {cmd}")
    subprocess.run(cmd, check=True, shell=True)


def hisat2(ubams, bam, index, rna_strandness, summary_file, threads, temp_dir):
    r1_reads = []
    r2_reads = []
    with tempfile.TemporaryDirectory(dir=temp_dir) as tmpdir:
        logger.debug(f"temp directory = {tmpdir}")
        for ubam in ubams:
            name = re.sub(".unmapped.bam$", "", os.path.basename(ubam))
            r1 = os.path.join(tmpdir, f"{name}_1.fq.gz")
            r2 = os.path.join(tmpdir, f"{name}_2.fq.gz")
            os.mkfifo(r1)
            r1_reads.append(r1)
            os.mkfifo(r2)
            r2_reads.append(r2)
            run(f"samtools fastq -1 {r1} -2 /dev/null -0 /dev/null -s /dev/null -O -n -F 0x900 {ubam} 2>/dev/null &")
            run(f"samtools fastq -1 /dev/null -2 {r2} -0 /dev/null -s /dev/null -O -n -F 0x900 {ubam} 2>/dev/null &")

        cmd = f"""hisat2 -x {index} -1 {",".join(r1_reads)} -2 {",".join(r2_reads)}"""
        if threads:
            cmd += f" -p {threads}"
        if rna_strandness:
            cmd += f" --rna-strandness {rna_strandness}"
        if summary_file:
            cmd += f" --summary-file {summary_file} --new-summary"

        cmd += f" --rg-id A,B,C "

        cmd += f" 2>/dev/null | samtools sort -o {bam}"
        logger.debug(cmd)
        run(cmd)
        cmd = f"samtools index -b {bam}"
        run(cmd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", dest="bam", required=True)
    parser.add_argument("-i", "--index", required=True)
    parser.add_argument("-s", "--rna-strandness", choices=["F", "R", "FR", "RF"])
    parser.add_argument("-S", "--summary-file")
    parser.add_argument("-p", dest="threads")
    parser.add_argument("-T", "--temp-dir", default=tempfile.gettempdir())
    parser.add_argument("ubams", nargs="+")
    args = parser.parse_args()

    ubams = []
    for ubam in args.ubams:
        if not ubam.endswith(".unmapped.bam"):
            logger.error(f"uBAM file name must end with .unmapped.bam was {ubam}")
            sys.exit(1)
        logger.info(f"Input uBAM: {ubam}")
        ubams.append(os.path.abspath(ubam))

    hisat2(
        ubams,
        os.path.abspath(args.bam),
        os.path.abspath(args.index),
        args.rna_strandness,
        args.summary_file,
        args.threads,
        args.temp_dir
    )


if __name__ == "__main__":
    main()
