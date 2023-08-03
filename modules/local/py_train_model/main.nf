// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRAIN_MODEL {

    label 'memory_medium'
    label 'gpu'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.4'

    input:
    path('assay_reshaped.h5')
    path('train.pb')
    path('test.pb')
    path('model.py')
    path('params.py')
    val(name)

    output:
    path("*history_model.pb"), emit: history
    path("*labels.pb"), emit: labels
    path(name), emit: model

    script:
    """
    train_model.py $name $task.cpus
    """
}
