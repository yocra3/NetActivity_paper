// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRAIN_TCGA {

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.4'

    input:
    path('assay_reshaped.h5')
    path('train.pb')
    path('test.pb')
    path('network_config.py')
    val(name)

    output:
    path("history*.pb"), emit: history
    path("model.pb"), emit: model

    script:
    """
    train_tcga.py network_config.py $name
    """
}
