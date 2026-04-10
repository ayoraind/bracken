// include { check_mandatory_parameter; check_optional_parameters; check_parameter_value } from './params_utilities.nf'

def default_params(){
    /***************** Setup inputs and channels ************************/
    def params = [:]
// Defaults for configurable variables
    params.help = false
    params.version = false
    params.database = false
    params.taxon = false
    params.sequencing_date = false
    params.kraken_report = false
    params.output_dir= "./.default_output"
    params.pipeline_info = "${params.output_dir}/pipeline_info"
    params.taxon = "S"
    
    return params
}

def check_mandatory_parameter(Map params, String parameter_name){
    if ( !params[parameter_name]){
        println "You must specify " + parameter_name
        System.exit(1)
    } else {
        return params[parameter_name]
    }
}

def check_optional_parameters(Map params, List parameter_names){
    if (parameter_names.collect{name -> params[name]}.every{param_value -> param_value == false}){
        println "You must specifiy at least one of these options: " + parameter_names.join(", ")
        System.exit(1)
    }
}

def check_parameter_value(String parameter_name, String value, List value_options){
    if (value_options.any{ it == value }){
        return value
    } else {
        println "The value for " + parameter_name + " must be one of " + value_options.join(", ")
        System.exit(1)
    }
}

def check_params(Map params) { 
    def final_params = params
    
    // set up fasta files
    final_params.assemblies = check_mandatory_parameter(params, 'assemblies') - ~/\/$/
     
     // set up output directory (use default if not provided)
    final_params.output_dir = params.containsKey('output_dir') ? params.output_dir : "./.default_output"
    final_params.output_dir = final_params.output_dir - ~/\/$/

    // set up output directory
    // final_params.output_dir = check_mandatory_parameter(params, 'output_dir') - ~/\/$/
     
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

