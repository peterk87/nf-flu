/*
 * This file holds several functions specific to the pipeline.
 */

class Workflow {

    // Citation string
    private static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
               "* The pipeline\n" + 
               "  TODO: add Zenodo DOI for pipeline\n\n" +
               "* The nf-core framework\n" +
               "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
               "* Software dependencies\n" +
               "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    static void validate_params(params, log, valid_params) {
        genome_exists(params, log)

        // Generic parameter validation
        if (!params.fasta) {
            log.error "Genome fasta file not specified!"
            System.exit(0)
        }
    }

    // Exit pipeline if incorrect --genome key provided
    static void genome_exists(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                      "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                      "  Currently, the available genome keys are:\n" +
                      "  ${params.genomes.keySet().join(", ")}\n" +
                      "==================================================================================="
            System.exit(0)
        }
    }

    // Print warning if genome fasta has more than one sequence
    static void is_multifasta(fasta, log) {
        def count = 0
        def line  = null
        fasta.withReader { reader ->
            while (line = reader.readLine()) {
                if (line.contains('>')) {
                    count++
                    if (count > 1) {
                        log.warn "=============================================================================\n" +
                                "  This pipeline does not officially support multi-fasta genome files!\n\n" + 
                                "  The parameters and processes are tailored for viral genome analysis.\n" +
                                "  Please amend the '--fasta' parameter.\n" +
                                "==================================================================================="
                        break
                    }
                }
            }
        }
    }

    // Function that parses and returns the number of mapped reasds from flagstat files
    static ArrayList get_flagstat_mapped_reads(workflow, params, log, flagstat) {
        def mapped_reads = 0
        flagstat.eachLine { line ->
            if (line.contains(' mapped (')) {
                mapped_reads = line.tokenize().first().toInteger()
            }
        }
        
        def pass = false
        def logname = flagstat.getBaseName() - 'flagstat'
        Map colors = Headers.log_colours(params.monochrome_logs)
        if (mapped_reads <= params.min_mapped_reads.toInteger()) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} [FAIL] Mapped read threshold >= ${params.min_mapped_reads}. IGNORING FOR FURTHER DOWNSTREAM ANALYSIS: ${mapped_reads} - $logname${colors.reset}."
        } else {
            pass = true
        }
        return [ mapped_reads, pass ]
    }
}
