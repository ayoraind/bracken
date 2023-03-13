params.kraken_report = "$projectDir/data/ggal/gut_{1,2}.fastq"
params.output_dir = "$PWD/results"
params.database = "/KRAKEN_DB/"
params.taxon = "S"
nextflow.enable.dsl=2




process BRACKEN_BRACKEN {
    tag "$meta"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI.
    // Please update version string below when bumping container versions.
    conda "bioconda::bracken=2.7"
    
    input:
    tuple val(meta), path(kraken_report)
    path database

    output:
    tuple val(meta), path(bracken_report), emit: reports
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta}"
    bracken_report = "${prefix}.tsv"
    // WARN: Version information not provided by tool on CLI.
    // Please update version string below when bumping container versions.
    def VERSION = '2.7'
    """
    bracken \\
        ${args} \\
        -d '${database}' \\
        -i '${kraken_report}' \\
        -o '${bracken_report}'
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: ${VERSION}
    END_VERSIONS
    """
}

workflow  {
         kraken_ch = channel
                          .fromPath( params.kraken_report, checkIfExists: true )
                          .map { file -> tuple(file.simpleName, file) }

        // taxon_ch = channel.of(params.taxon)

        // BRACKEN(kraken_ch, params.database, taxon_ch)
        // BRACKEN(kraken_ch, params.database)
	 BRACKEN_BRACKEN(kraken_ch, params.database)
}
