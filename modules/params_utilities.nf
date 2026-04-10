include { version_message; help_message } from './messages.nf'

def help_or_version(Map params, String version){
    // Show help message
    if (params.help){
        version_message(version)
        help_message()
        System.exit(0)
    }

    // Show version number
    if (params.version){
        version_message(version)
        System.exit(0)
    }
}


def rename_params_keys(Map params_to_rename, Map old_and_new_names) {
    old_and_new_names.each{ old_name, new_name ->
        if (params_to_rename.containsKey(old_name))  {
            params_to_rename.put( new_name, params_to_rename.remove(old_name ) )
        }
    }
    return params_to_rename
}
