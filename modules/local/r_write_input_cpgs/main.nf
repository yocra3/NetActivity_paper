// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process R_WRITE_INPUT_CPGS {

    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_rsession:1.1'

    input:
    tuple val(prefix), path(hdf5), path(rds)

    output:
    path "*.txt", emit: res

    script:
    """
    write_input_CpGs.R $hdf5
    """
}
