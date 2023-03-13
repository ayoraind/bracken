nextflow.enable.dsl=2

params.kraken_report = "$projectDir/data/ggal/gut_{1,2}.fastq"
params.reads = "$projectDir/data/ggal/gut_{1,2}.fastq"
params.output_dir = "$PWD/results"
nextflow.enable.dsl=2
params.database = "/KRAKEN_DB/"
params.taxon = "S"


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

  tag { "abundance estimation in ${sample_id}" }
    
  if (params.output_dir != "") { publishDir(params.output_dir, mode:'copy') }
  
  conda "bioconda::bracken=2.8"
  
  input:
  tuple val(sample_id), path(kraken_report)
  path(brackendb)
  val(taxon)
  

  output:
  path("*_bracken"),                   emit: bracken_out_ch
  path("*_bracken_kraken.report.txt"), emit: bracken_kraken_report_ch
  path("*bracken.log.file.txt"),       emit: bracken_log_ch
  path "versions.yml", 		       emit: versions
 
  
  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${sample_id}"
  def VERSION = '2.8'
  """
      
  bracken ${args} -d ${brackendb} -i ${kraken_report} -o ${prefix}.${taxon}_bracken -r 100 -l ${taxon} -t ${task.cpus} -w ${prefix}_${taxon}_bracken_kraken.report.txt > ${prefix}.bracken.log.file.txt
  
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
          //kraken_ch.view()
	  
	  database_ch = channel.fromPath( params.database ) 
	  //database_ch.view()
	  
	  BRACKEN(kraken_ch, database_ch.collect(), params.taxon)
       
         //BRACKEN(kraken_ch, params.database, params.taxon)
}
