// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process EXTRACT_MODEL_FEATURES {

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.4'

    input:
    path('assay_reshaped.h5')
    path('model/')

    output:
    path("*.tsv"), emit: res

    script:
    """
    extract_features.py
    """
}
