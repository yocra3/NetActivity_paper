// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRAIN_AUTOENCODER_MODEL {

    label 'memory_high'
    label 'cpus_high'
    label 'process_long'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.3'

    input:
    path('train.pb')
    path('test.pb')
    path('model.py')
    path('params.py')
    // Old code for autoencoder after CNN
    //path('model/')
    val(name)

    output:
    // path("*history_model.pb"), emit: old_history
    // path("*autoencoder_history_model.pb"), emit: auto_history
    // path("*_with_autoencode"), emit: model_old
    // path("*_autoencoder"), emit: model_auto
    path("*history_model.pb"), emit: history
    path("*_autoencoder"), emit: model

    script:
    """
    train_autoencoder.py $name $task.cpus
    """
}
