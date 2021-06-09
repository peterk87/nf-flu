/*
 * Helper functions for peterk87/nf-virontus
 */


//=============================================================================
// User input validation helper functions
//=============================================================================

def checkFileExists(file_path) {
  f = file(file_path)
  if ( !f.isFile() || !f.exists() ) {
    exit 1, "File '$file_path' does not exist!"
  }
}

def checkCentrifugeDb(centrifuge_db) {
  file_centrifuge_db = file(centrifuge_db)
  prefix = file_centrifuge_db.getName()
  centrifuge_dir = file_centrifuge_db.getParent()
  if ( !centrifuge_dir.isDirectory() || !centrifuge_dir.exists() ) {
    exit 1, "Centrifuge DB does not exist at '$centrifuge_dir'! Please specify a valid Centrifuge DB."
  }
  any_valid = false
  centrifuge_dir.eachFile { f ->
    if ( f.isFile() ) {
      if ( f.getName() =~ /^$prefix/ && f.getExtension() == 'cf') {
        any_valid = true
      }
    }
  }
  if ( !any_valid ) {
    exit 1, "No valid Centrifuge DB files with prefix '$prefix' in '$centrifuge_dir' and extension 'cf'! Please specify a valid Centrifuge classification DB directory and prefix."
  }
}

def checkKraken2Db(kraken2_db) {
  kraken2_db_dir = file(kraken2_db)
  if ( !kraken2_db_dir.isDirectory() ) {
    exit 1, "The Kraken2 DB must be a directory! '$kraken2_db' is not a directory!"
  }
  if ( !kraken2_db_dir.exists() ) {
    exit 1, "The Kraken2 DB must be an existing directory! '$kraken2_db' does not exist!"
  }
}

// Check that all taxids are integers delimited by commas
def checkTaxids(taxids) {
  if (taxids instanceof Boolean || taxids.toString().isEmpty()) {
    return null
  } else if (taxids.toString().isInteger()) {
    return taxids.toString()
  } else {
    taxids_list = taxids.toString()
      .split(',')
      .collect { it.strip() }
      .findAll { it != '' }
    if (!taxids_list.every { it.isInteger() }) {
      exit 1, "Not every element in `--taxids` is an integer!"
    }
    return taxids_list.join(',')
  }
}

def check_sample_sheet(LinkedHashMap sample_sheet) {
  // Check that each entry from a sample sheet
  reads_file = file(sample_sheet.reads, checkIfExists: true)
  reads = sample_sheet.reads ? reads_file : null
  if (reads == null) {
    exit 1, "The Nanopore reads FASTQ file or directory specified for ${sample_sheet.sample} does not exist! Please check the sample sheet '$params.sample_sheet'"
  }
  if (reads_file.isDirectory()) {
    fqs = []
    reads_file.eachFile {
      fname = it.getName()
      if (fname ==~ /.*\.fastq(\.gz)?$/) {
        fqs << it
      }
    }
    if (fqs.size() == 0) {
      exit 1, "Sample '${sample_sheet.sample}' input specified as directory '$reads_file' with no FASTQ files! Please check the sample sheet '$params.sample_sheet'\n${fqs}\n${reads_file.listFiles()}"
    }
    reads = fqs
  }
  return [ sample_sheet.sample, reads ]
}
