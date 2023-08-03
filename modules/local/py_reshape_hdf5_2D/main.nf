// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PY_RESHAPE_HDF5_2D {

    label 'process_high'
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
    reshape_hdf5_2D.py
    """
}
