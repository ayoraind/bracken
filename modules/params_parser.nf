include { check_mandatory_parameter; check_optional_parameters; check_parameter_value } from './params_utilities.nf'

def default_params(){
    /***************** Setup inputs and channels ************************/
    def params = [:] as nextflow.script.ScriptBinding$ParamsMap
    // Defaults for configurable variables
    params.help = false
    params.version = false
    params.database = false
    params.taxon = false
    params.sequencing_date = false
    params.kraken_report = false
    return params
}

def check_params(Map params) { 
    final_params = params
     
    // set up output directory
    final_params.output_dir = check_mandatory_parameter(params, 'output_dir') - ~/\/$/
     
    // set up database directory
    final_params.database = check_mandatory_parameter(params, 'database') - ~/\/$/
    
    // set up sequencing date
    final_params.sequencing_date = check_mandatory_parameter(params, 'sequencing_date')
    
    // set up kraken report
    final_params.kraken_report = check_mandatory_parameter(params, 'kraken_report')
	
    // set up taxon
    final_params.taxon = check_mandatory_parameter(params, 'taxon') 
        
    // check taxon is valid
    final_params.taxon = check_parameter_value('taxon', final_params.taxon, ['S', 'S1', 'G', 'G1', 'F', 'O', 'C', 'P', 'D'])
    
    return final_params
}

