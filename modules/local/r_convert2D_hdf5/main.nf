// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process R_CONVERT2D_HDF5 {

    label 'process_high'
    label 'process_long'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_rsession:1.1'

    input:
    tuple val(prefix), path(hdf5), path(rds)

    output:
    path "*.h5", emit: res

    script:
    """
    convert2D.R $prefix
    """
}
