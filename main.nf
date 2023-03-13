#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include definitions
include  { helpMessage; Version } from './modules/messages.nf'

// include processes
include { BRACKEN; EXTRACT_TAXON_SPECIFIC_INFO; COMBINE_BRACKEN_KRAKEN_REPORT_FROM_TAXA; COMBINE_BRACKEN_FILES_FROM_TAXA;  EXTRACT_UNCLASSIFIED_INFO_FROM_BRACKEN_LOG; COMBINE_UNCLASSIFIED_LOGS } from './modules/processes.nf'

log.info """\
    ======================================
    BRACKEN  - TAPIR  P I P E L I N E
    ======================================
    output_dir      : ${params.output_dir}
    """
    .stripIndent()


workflow  {
          kraken_ch = channel
                          .fromPath( params.kraken_report, checkIfExists: true )
                          .map { file -> tuple(file.simpleName, file) }


          database_ch = channel.fromPath( params.database )


          BRACKEN(kraken_ch, database_ch.collect(), params.taxon)

          //for bracken kraken reports
          EXTRACT_TAXON_SPECIFIC_INFO(BRACKEN.out.bracken_kraken_report_ch, params.taxon)

          collected_bracken_kraken_taxon_ch = EXTRACT_TAXON_SPECIFIC_INFO.out.taxon_bracken_kraken_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )

          COMBINE_BRACKEN_KRAKEN_REPORT_FROM_TAXA(collected_bracken_kraken_taxon_ch, params.taxon, params.sequencing_date)

          //for bracken reports
          collected_bracken_taxon_ch = BRACKEN.out.bracken_out_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )

          COMBINE_BRACKEN_FILES_FROM_TAXA(collected_bracken_taxon_ch, params.taxon, params.sequencing_date)

          //for extracting unclassified read info from bracken log

          EXTRACT_UNCLASSIFIED_INFO_FROM_BRACKEN_LOG(BRACKEN.out.bracken_log_ch)

          collected_unclassified_log_ch = EXTRACT_UNCLASSIFIED_INFO_FROM_BRACKEN_LOG.out.bracken_uncl_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )
          COMBINE_UNCLASSIFIED_LOGS(collected_unclassified_log_ch, params.sequencing_date)
}
