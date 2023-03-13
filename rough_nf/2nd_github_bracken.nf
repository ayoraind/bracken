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
  tuple val(sample_id), path(reads)
  val(brackendb)
  file(kraken_report)
  val(taxon)
  

  output:
  path("*_bracken"),                   emit: bracken_out_ch
  path("*_bracken_kraken.report.txt"), emit: bracken_kraken_report_ch
  path("*_bracken_log_file.txt"),      emit: bracken_log_ch
  path "versions.yml", 		       emit: versions
 
  
  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${sample_id}"
  def VERSION = '2.8'
  """
  TAXON_UPPERCASE=\$(echo ${taxon} | tr '[:lower:]' '[:upper:]')
  FIRST=\$(echo ${reads} | head -n1 | awk '{print \$1;}')
  if [ `echo \$FIRST | grep "gz"` ]; then cat="zcat"; else cat="cat"; fi
  READSIZE=\$(\$cat \$FIRST | awk '{num++}{if (num%4==2){line++; sum+=length(\$0)} if (line==100) {printf "%.0f", sum/100; exit} } ')
    
  bracken ${args} -d ${brackendb} -i ${kraken_report} -o ${prefix}.\${TAXON_UPPERCASE}_bracken -r \$READSIZE -l \${TAXON_UPPERCASE} -t ${task.cpus} -w ${prefix}_\${TAXON_UPPERCASE}_bracken_kraken.report.txt > ${prefix}.bracken.log.file.txt
  
  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: ${VERSION}
  END_VERSIONS
   
  """

  }


workflow  {
         reads_ch = channel
                          .fromPath( params.reads, checkIfExists: true )
                          .map { file -> tuple(file.simpleName, file) }

       
         BRACKEN(reads_ch, params.database, params.kraken_report, params.taxon)
}
