module main

import time
import os
import flag
import strings
import encoding.csv

const segments = {
	'A_PB2': '1_PB2'
	'A_PB1': '2_PB1'
	'A_PA':  '3_PA'
	'A_HA':  '4_HA'
	'A_NP':  '5_NP'
	'A_NA':  '6_NA'
	'A_MP':   '7_M'
	'A_NS':  '8_NS'
	'B_PB1': '1_PB1'
	'B_PB2': '2_PB2'
	'B_PA':  '3_PA'
	'B_HA':  '4_HA'
	'B_NP':  '5_NP'
	'B_NA':  '6_NA'
	'B_MP':   '7_M'
	'B_NS':  '8_NS'
}

fn refname_to_segment(x string) string {
	mut segment_key := x
	if x.starts_with('A_HA') || x.starts_with('A_NA') {
		segment_parts := x.split('_')
		if segment_parts.len >= 2 {
			segment_prefix := segment_parts[0] + '_' // "A_"
			segment_suffix := segment_parts[1] // "HA", "NA", etc.
			segment_key = segment_prefix + segment_suffix // Reconstructs the segment identifier
		}
	}
	return segments[segment_key] or { x }
}

fn parse_tsv_to_fasta(input_file string, output_file string) ! {
	mut consensus_sequence := map[int]string{}
	mut position_counts := map[int]map[string]int{}

	mut reader := csv.csv_sequential_reader(separator: `\t`, file_path: input_file)!
	header := reader.get_next_row()!
	pos_idx := header.index('Position')
	allele_idx := header.index('Allele')
	count_idx := header.index('Count')
	for {
		row := reader.get_next_row() or { break }
		if row.len == 0 {
			break
		}
		position := row[pos_idx].int()
		allele := row[allele_idx]
		if allele == '-' {
			continue
		}
		count := row[count_idx].int()

		if position_counts[position].len == 0 {
			position_counts[position] = map[string]int{}
		}
		position_counts[position][allele] = count
	}

	for position, alleles in position_counts {
		mut max_count := 0
		mut consensus_allele := ''
		for allele, count in alleles {
			if count > max_count {
				consensus_allele = allele
				max_count = count
			}
		}
		consensus_sequence[position] = consensus_allele
	}

	mut sorted_positions := consensus_sequence.keys()
	sorted_positions.sort()
	mut sequence := strings.new_builder(sorted_positions.len)
	for pos in sorted_positions {
		sequence.write_string(consensus_sequence[pos])
	}

	mut outfile := os.create(output_file) or { panic(err) }
	defer {
		outfile.close()
	}
	base_name := os.base(output_file)
	name := base_name.all_before_last('.')
	outfile.writeln('>${name}')!
	outfile.writeln(sequence.str())!
}

fn main() {
	mut fp := flag.new_flag_parser(os.args)
	fp.application('Parse TSV files to FASTA format.')
	fp.skip_executable()

	sample_name := fp.string('sample-name', `n`, '', 'Sample name for output file naming.')
	input_directory := fp.string('input-directory', `i`, '', 'Directory containing input TSV files.')
	output_directory := fp.string('output-directory', `o`, '', 'Directory to save output FASTA files.')
	fp.finalize() or { panic(err) }
	if sample_name == '' {
		eprintln('Sample name needs to be specified')
		exit(1)
	}
	if input_directory == '' {
		eprintln('Input directory needs to be specified')
		exit(1)
	}
	if !os.exists(input_directory) {
		eprintln('Input directory ${input_directory} does not exist!')
		exit(1)
	}
	if output_directory == '' {
		eprintln('Output directory needs to be specified')
		exit(1)
	}

	if !os.exists(output_directory) {
		os.mkdir_all(output_directory) or { panic(err) }
	}

	sw := time.new_stopwatch()
	files := os.ls(input_directory) or { panic(err) }
	for file in files {
		if file.ends_with('-allAlleles.txt') {
			eprintln('Parsing ${file}')
			input_file := os.join_path(input_directory, file)
			segment := refname_to_segment(file.all_before_last('-allAlleles'))
			output_file := os.join_path(output_directory, '${sample_name}_${segment}.fasta')
			parse_tsv_to_fasta(input_file, output_file) or { panic(err) }
		}
	}
	eprintln('Done in ${sw.elapsed().milliseconds()} ms')
}
