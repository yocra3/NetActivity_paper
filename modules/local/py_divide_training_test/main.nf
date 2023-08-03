// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PY_DIVIDE_TRAIN_TEST {

    label 'memory_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.4'

    input:
    path('assays.h5')
    val(prop)
    val(autoencoder)

    output:
    path("train.pb"), emit: train
    path("test.pb"), emit: test
    path("test_indices.csv"), emit: indices
    
    script:
    """
    divide_train_test.py $prop $autoencoder
    """
}
