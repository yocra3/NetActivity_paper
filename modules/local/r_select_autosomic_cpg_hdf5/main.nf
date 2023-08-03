// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process R_SELECT_AUTOSOMICCPG_HDF5 {

    label 'memory_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_rsession:1.1'

    input:
    tuple val(prefix), path(hdf5), path(rds)

    output:
    tuple val("${prefix}autosomicProbes_"), path("*.h5"), path("*.rds"), emit: res

    script:
    """
    select_autosomicCpG_hdf5.R $hdf5 $prefix
    """
}
