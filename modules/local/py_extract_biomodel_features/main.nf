// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PY_EXTRACT_BIOMODEL_FEATURES {

    label 'memory_medium'
    label 'process_long'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.3b'

    input:
    path('assays.h5')
    path('gene_mask.pb')
    path(checkpoints)
    path('model.py')
    path('params.py')

    output:
    path("*.tsv"), emit: res

    script:
    """
    extract_features_bio.py
    """
}
