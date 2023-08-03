// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PY_TRAIN_DNN_BIOLOGICALMODEL {

    label 'memory_medium'
    label 'cpus_medium'
    label 'process_long'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.3b'

    input:
    path('assay_reshaped.h5')
    path('train.pb')
    path('test.pb')
    path('model.py')
    path('params.py')
    path('gene_mask.pb')
    val(name)


    output:
    path("*history_model.pb"), emit: history
    path("*labels.pb"), emit: labels
    path('checkpoints'), emit: checkpoints


    script:
    """
    train_dnn_biologicalmodel.py $name $task.cpus
    """
}
