#!/usr/bin/env python
"""
This script is modified from nf-core/viralrecon script:
https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/bin/fastq_dir_to_samplesheet.py#L1

Modifications:
- regex to extract the sample name from FASTQ filenames instead of stripping away forward/reverse read file suffixes.
- remove "Undetermined_*.fastq.gz" sample and reads from samplesheet
- defaults for input dir (current directory: '.') and output samplesheet ("samplesheet.csv")
"""

import argparse
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path


def parse_args(args=None) -> argparse.Namespace:
    desc = "Generate a samplesheet.csv from a directory of FastQ files."
    epilog = "Example usage: python fastq_dir_to_samplesheet.py -i <FASTQ_DIR> -o <SAMPLESHEET_FILE>"

    parser = argparse.ArgumentParser(description=desc, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input-dir",
        type=str,
        default=".",
        help="Directory containing Illumina FASTQ files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="samplesheet.csv",
        help="Output samplesheet CSV file.",
    )
    parser.add_argument(
        "--sample-name-regex",
        type=str,
        default=r"(.+)_S\d+_L\d{3}_R[12]_001\.fastq\.gz",
        help="Sample name regular expression where first matching group is the desired sample name.",
    )
    parser.add_argument(
        "-f", "--force", action="store_true", help="Overwrite samplesheet CSV?"
    )
    parser.add_argument(
        "--keep-undetermined",
        action="store_true",
        help='Keep the "Undetermined*R1/2*.fastq.gz" files? Default is to remove these reads from the samplesheet.',
    )
    return parser.parse_args(args)


def fastq_dir_to_samplesheet(fastq_dir: Path,
                             samplesheet_file: Path,
                             sample_name_regex: str,
                             keep_undetermined: bool = False) -> None:
    sample_name_regex = re.compile(sample_name_regex)
    read_dict = defaultdict(list)
    for path in fastq_dir.glob("*"):
        m = sample_name_regex.match(path.name)
        if m:
            sample = m.group(1)
            read_dict[sample].append(str(path.absolute()))
    if not keep_undetermined and 'Undetermined' in read_dict:
        logging.info(f'Removing Undetermined FASTQ files from samplesheet.')
        del read_dict['Undetermined']

    logging.info(f'Found {len(read_dict)} samples with {sum(len(v) for v in read_dict.values())} FASTQ files.')
    if read_dict:
        with open(samplesheet_file, "w") as fout:
            header = ["sample", "fastq_1", "fastq_2"]
            fout.write(",".join(header) + "\n")
            for sample, reads in read_dict.items():
                sample_info = ",".join([sample] + reads)
                if len(reads) == 1:
                    sample_info += ","
                fout.write(f"{sample_info}\n")
        logging.info(f'Wrote samplesheet to "{samplesheet_file}".')
    else:
        error_str = (
            "No FastQ files found so samplesheet has not been created!"
            "Please check the values provided for the:\n"
            "  - Path to the directory containing the FastQ files\n"
            "  - '--read1_extension' parameter\n"
            "  - '--read2_extension' parameter"
        )
        logging.error(error_str)
        sys.exit(1)


def main(args=None):
    args = parse_args(args)

    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]",
        level=logging.INFO,
    )

    inputdir = Path(args.input_dir)
    if not (inputdir.exists() and inputdir.is_dir()):
        raise NotADirectoryError(
            f'Input FASTQ directory "{inputdir}" does not exist! '
            f"Please provide an input directory that exists to `-i/--inputdir`"
        )

    output = Path(args.output)
    if output.exists():
        if not args.force:
            raise FileExistsError(
                f'Samplesheet file already exists at "{output}"! Overwrite with `-f/--force`'
            )
        else:
            logging.warning(f'Overwriting existing samplesheet at "{output}"!')

    fastq_dir_to_samplesheet(
        fastq_dir=inputdir,
        samplesheet_file=output,
        sample_name_regex=args.sample_name_regex,
    )


if __name__ == "__main__":
    sys.exit(main())
