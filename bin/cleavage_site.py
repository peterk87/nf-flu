#!/usr/bin/env python

import logging
import click
import re
from Bio import SeqIO

LOG_FORMAT = "%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]"
logging.basicConfig(format=LOG_FORMAT, level=logging.INFO)

def extract_sequences(input_file):
    ha1_ha2_pairs = []

    ha1_sequence = None
    ha2_sequence = None
    common_id = None

    for seq_record in SeqIO.parse(input_file, "fasta"):
        full_header = seq_record.description.replace(" ", "_")
        sequence = str(seq_record.seq)

        if "segment4_HA" in full_header:
            logging.info(f"Full Header: {full_header}")
            logging.info(f"First 30 characters of Sequence: {sequence[:30]}...")

            if "HA1" in full_header:
                ha1_sequence = sequence[-15:]
                common_id = full_header.split("_")[0]  # Extract the common identifier from the header
                logging.info(f"HA1 sequence extracted: {ha1_sequence}")

            elif "HA2" in full_header:
                if common_id and full_header.split("_")[0] == common_id:
                    ha2_sequence = sequence[:3]
                    logging.info(f"HA2 sequence extracted: {ha2_sequence}")

                    if ha1_sequence and ha2_sequence:
                        ha1_ha2_pairs.append((ha1_sequence, ha2_sequence))
                        logging.info(f"Recognized HA1/HA2 pair: {ha1_sequence}|{ha2_sequence}")

    logging.info(f"Extracted {len(ha1_ha2_pairs)} HA1/HA2 pairs")

    return ha1_ha2_pairs

def process_cleavage_site(output_file, ha1_ha2_pairs):
    """Process the sequences to detect cleavage sites and classify them."""
    logging.info(f"Writing results to: {output_file}")

    # Log the HA1/HA2 pairs received by process_cleavage_site
    logging.info(f"HA1/HA2 pairs received in process_cleavage_site: {ha1_ha2_pairs}")

    # Open the output file for writing results
    with open(output_file, 'w') as output:
        output.write("Cleavage Sequence\tResidue Type\tBasic Residue Count\tPathogenicity\n")

        # If no pairs were passed, log and exit
        if not ha1_ha2_pairs:
            logging.error("No cleavage sequences found. Please check the input file or extraction logic.")
            return

        for idx, (ha1_sequence, ha2_sequence) in enumerate(ha1_ha2_pairs):
            logging.info(f"Processing pair {idx + 1}:")
            logging.info(f"HA1 sequence: {ha1_sequence}")
            logging.info(f"HA2 sequence: {ha2_sequence}")

            # Regular expression for cleavage site detection in HA1 sequence
            ha1_pattern = r"P(Q|L|E)[A-Z]{1,10}(R|K)(.*)"
            ha1_match = re.search(ha1_pattern, ha1_sequence)
            if ha1_match:
                ha1_sequence = ha1_match.group(0)
                logging.info(f"Cleavage site detected in HA1: {ha1_sequence}")
            else:
                ha1_sequence = "No cleavage site found"
                logging.info("No cleavage site detected in HA1")

            # Combine HA1 and HA2 to form the cleavage sequence
            cleavage_sequence = f"{ha1_sequence}|{ha2_sequence}"
            logging.info(f"Cleavage sequence formed: {cleavage_sequence}")

            # Remove 'P.(R|K)' pattern from HA1 for residue counting
            ha1_sequence = re.sub(r"P.(R|K)", "", ha1_sequence)

            # Initialize variables for residue counting
            basic_residue_count = 0
            stop_counting = False
            residue_type = "No match"

            # Count basic residues in HA1 sequence
            if ha1_sequence[-1] not in "RK":
                logging.info("Position -1 is not basic! Labeling as 'check seq'.")
                residue_type = "check seq"
                basic_residue_count = 0
                stop_counting = True
            else:
                consecutive_non_basic = 0
                for i in range(-1, -len(ha1_sequence) - 1, -1):
                    if ha1_sequence[i] in "RK":
                        basic_residue_count += 1
                        consecutive_non_basic = 0
                    else:
                        consecutive_non_basic += 1
                        if abs(i) == 2 or abs(i) == 3:
                            if consecutive_non_basic > 1:
                                stop_counting = True
                                break
                            continue
                        else:
                            stop_counting = True
                            break

                # Classify based on basic residue count
                if stop_counting or i == -len(ha1_sequence):
                    if basic_residue_count == 1:
                        residue_type = "Monobasic"
                    elif basic_residue_count > 1:
                        residue_type = "Multibasic"
                    else:
                        residue_type = "No match"

            # Determine pathogenicity
            if basic_residue_count == 1:
                pathogenicity = "LPAI"
            elif basic_residue_count == 2:
                pathogenicity = "LP/HP"
            elif basic_residue_count >= 3:
                pathogenicity = "HPAI"
            else:
                pathogenicity = "check seq"

            # Write the results to the output file in TSV format
            output.write(f"{cleavage_sequence}\t{residue_type}\t{basic_residue_count}\t{pathogenicity}\n")
            logging.info(f"Cleavage sequence output for pair {idx + 1}: {cleavage_sequence}, Type: {residue_type}, Basic residue count: {basic_residue_count}, Pathogenicity: {pathogenicity}")


# Main script
@click.command()
@click.option("-i", "--input-file", default="input.faa", help="Input FASTA file")
@click.option("-o", "--output-file", default="ha_cleavage_site.tsv", help="Output file for cleavage site analysis")
@click.option("-s", "--sample-name", default=None, help="Sample name to use")
def main(input_file, output_file, sample_name):
    """Main function to process the FASTA file and extract cleavage sequences."""
    # Step 1: Extract HA1 and HA2 sequences
    ha1_ha2_pairs = extract_sequences(input_file)

    # Step 2: Process cleavage site analysis and classification
    process_cleavage_site(output_file, ha1_ha2_pairs)


if __name__ == "__main__":
    main()
