// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PY_COMBINE_METHY_LABEL_HDF5 {

    label 'memory_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.4'

    input:
    path('assays.h5')
    path('labels.txt')

    output:
    path("*.h5"), emit: assay

    script:
    """
    combine_methy_label_hdf5.py
    """
}
