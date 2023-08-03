// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRANSFER_LEARNING_DENSE {

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.4'

    input:
    path('assay_reshaped.h5')
    path('train.pb')
    path('test.pb')
    path('model.pb')
    val(name)

    output:
    path("history*.pb"), emit: history
    path("new_model.pb"), emit: model

    script:
    """
    transfer_learning_dense.py network_config.py $name
    """
}
