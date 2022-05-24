//Credit https://github.com/nf-core/viralrecon/blob/master/modules/local/multiqc_tsv_from_list.nf
process MULTIQC_TSV_FROM_LIST {

    executor 'local'
    memory 100.MB

    input:
    val tsv_data   // [ ['foo', 1], ['bar', 1] ]
    val header     // [ 'name', 'number' ]
    val out_prefix

    output:
    path "*.tsv"

    when:
    task.ext.when == null || task.ext.when

    exec:
    // Generate file contents
    def contents = ""
    if (tsv_data.size() > 0) {
        contents += "${header.join('\t')}\n"
        contents += tsv_data.join('\n')
    }

    // Write to file
    def mqc_file = task.workDir.resolve("${out_prefix}_mqc.tsv")
    mqc_file.text = contents
}
