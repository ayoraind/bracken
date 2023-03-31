/*
 * pipeline input parameters
 */

params.kraken_report = "$projectDir/data/ggal/gut_{1,2}.fastq"
params.output_dir = "$PWD/results"
params.database = "/KRAKEN_DB/"
params.taxon = "S"
nextflow.enable.dsl=2


log.info """\
    BRACKEN  - TAPIR   P I P E L I N E
    ============================================
    output_dir       : ${params.output_dir}
    database         : ${params.database}
    taxon            : ${params.taxon}
    """
    .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
 
process BRACKEN {
    publishDir "${params.output_dir}", mode:'copy'
    tag "species abundance estimation from ${meta}"
    
    
    conda "bioconda::bracken=2.8"
    
    input:
    tuple val(meta), path(kraken_report)
    path(database)


    output:
    tuple val(meta), path(bracken_report),        emit: bracken_report_ch
    tuple val(meta), path(bracken_kraken_report), emit: bracken_kraken_report_ch
    path "versions.yml"                         , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
   // bracken_report = ${prefix}.${taxon}.bracken
    bracken_report = ${prefix}.bracken
    // bracken_kraken_report = ${prefix}.${taxon}.bracken.kraken.report.txt
    bracken_kraken_report = ${prefix}.bracken.kraken.report.txt
    // WARN: Version information not provided by tool on CLI.
    // Please update version string below when bumping container versions.
    def VERSION = '2.8'
    """
    # species abundance estimation
    bracken ${args} -d '${database}' -i '${kraken_report}' -o '${bracken_report}' -r 100 -l ${taxon} -w '${bracken_kraken_report}'
    # bracken ${args} -d '${database}' -i '${kraken_report}' -o ${prefix}.${taxon}.bracken -r 100 -l ${taxon} -w ${prefix}.${taxon}.bracken.kraken.report.txt
    
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
	 BRACKEN(kraken_ch, params.database)  
}

