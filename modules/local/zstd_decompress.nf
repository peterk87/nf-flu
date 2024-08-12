process ZSTD_DECOMPRESS {

  conda 'conda-forge::zstd=1.5.2'
  // TODO: using clair3 container here for zstd and since it might be used if running the Nanopore workflow, but should move to multi-package-container with just zstd and maybe curl to combine data fetch functionality
  container 'hkubal/clair3:v1.0.10'

  input:
  path(zstd_file, stageAs: "input*/*")
  val(filename)

  output:
  path(decompressed_file), emit: file
  path('versions.yml'), emit: versions

  script:
  def basename = file(zstd_file).getBaseName()
  decompressed_file = filename ? "${basename}-${filename}" : basename
  """
  zstdcat $zstd_file > $decompressed_file

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      zstd: \$(echo \$(zstd --version 2>&1) | sed 's/^.* v//; s/,.*//')
  END_VERSIONS
  """
}
