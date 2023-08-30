def help_message() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --reads "PathToReadFile(s)" --output_dir "PathToOutputDir" --database "PathToDB" --sequencing_date "GYYMMDD" --taxon "taxon_symbol" --kraken_report "*.kraken.report.txt"

        Mandatory arguments:
         --database                	KRAKEN database directory (full path required, e.g., "/KRAKEN_DB")
         --output_dir                   Output directory to place final combined bracken kraken output (e.g., "/MIGE/01_DATA/11_QC_REPORTS")
         --sequencing_date              Sequencing Date (for TAPIR, must start with G e.g., "G230223")
	 --kraken_report                kraken report file generated from kraken analysis

        Optional arguments:
         --taxon                        value for taxon must be one of: "D", "P", "C", "O", "F", "G", "G1", "S", "S1" (Default: "S")
         --help                         This usage statement.
         --version                      Version statement
        """
}



def version_message(String version) {
      println(
            """
            ===========================================
             BRACKEN TAPIR Pipeline version ${version}
            ===========================================
            """.stripIndent()
        )

}


def pipeline_start_message(String version, Map params){
    log.info "======================================================================"
    log.info "                  BRACKEN TAPIR Pipeline version ${version}           "
    log.info "======================================================================"
    log.info "Running version   : ${version}"
    log.info ""
    log.info "-------------------------- Other parameters --------------------------"
    params.sort{ it.key }.each{ k, v ->
        if (v){
            log.info "${k}: ${v}"
        }
    }
    log.info "======================================================================"
    log.info "Outputs written to path '${params.output_dir}'"
    log.info "======================================================================"

    log.info ""
}


def complete_message(Map params, nextflow.script.WorkflowMetadata workflow, String version){
    // Display complete message
    log.info ""
    log.info "Ran the workflow: ${workflow.scriptName} ${version}"
    log.info "Command line    : ${workflow.commandLine}"
    log.info "Completed at    : ${workflow.complete}"
    log.info "Duration        : ${workflow.duration}"
    log.info "Success         : ${workflow.success}"
    log.info "Work directory  : ${workflow.workDir}"
    log.info "Exit status     : ${workflow.exitStatus}"
    log.info ""
}

def error_message(nextflow.script.WorkflowMetadata workflow){
    // Display error message
    log.info ""
    log.info "Workflow execution stopped with the following message:"
    log.info "  " + workflow.errorMessage
}
