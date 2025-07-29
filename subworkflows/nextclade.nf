nextflow.enable.dsl = 2

include { NEXTCLADE_DATASET_GET; NEXTCLADE_RUN; AGG_NEXTCLADE_TSV } from '../modules/local/nextclade'

workflow NEXTCLADE {
  take:
  ch_input_fasta            // channel: [sample, path]
  nextclade_datasets_csv    // params.nextclade_dataset_csv

  main:
  ch_versions = Channel.empty()
  ch_nextclade_datasets = Channel.fromPath(nextclade_datasets_csv)
    .splitCsv(header: false, strip: true)
    .map { row ->
      if (row.size() == 0){
        return null
      }
      def dataset_name = row[0]
      def dataset_tag = row.size() > 1 ? row[1] : null
      return [name: dataset_name, tag: dataset_tag]
    }
    .filter { it != null }
  NEXTCLADE_DATASET_GET(ch_nextclade_datasets)
  ch_versions = ch_versions.mix(NEXTCLADE_DATASET_GET.out.versions)
  ch_datasets = NEXTCLADE_DATASET_GET.out.dir
    .map { dataset, dataset_dir ->
      // get the actual version tag from the pathogen.json file
      def jsonFile = file("${dataset_dir}/pathogen.json")
      def json = new groovy.json.JsonSlurper().parseText(jsonFile.text)
      dataset.tag = json.version.tag
      return [dataset, dataset_dir]
    }
  NEXTCLADE_RUN(ch_input_fasta.combine(ch_datasets))
  ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions)
  ch_nextclade_outputs_csv = NEXTCLADE_RUN.out.tsv
    .map { sample, dataset, nextclade_tsv ->
      def filename = nextclade_tsv.getName()
      return "${sample},${dataset.name},${dataset.tag},${filename}"
    }
    .collectFile(name: "nextclade-tsv-outputs.csv", newLine: true)
  AGG_NEXTCLADE_TSV(
    NEXTCLADE_RUN.out.tsv.map { it[2] }.collect(),
    ch_nextclade_outputs_csv
  )
  ch_versions = ch_versions.mix(AGG_NEXTCLADE_TSV.out.versions)
  
  emit:
  versions = ch_versions
}
