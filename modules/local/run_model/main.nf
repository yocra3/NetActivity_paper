// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RUN_MODEL {

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.4'

    input:
    path('model.pb')
    path('assay_reshaped.h5')


    output:
    path("*.tsv"), emit: res

    script:
    """
    run_model.py
    """
}
