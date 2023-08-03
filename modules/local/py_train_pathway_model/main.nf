// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRAIN_PATHWAY_MODEL {

    label 'memory_medium'
    label 'gpu'
    label 'process_long'
    label 'cpu_high'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.4'

    input:
    path('assay_reshaped.h5')
    path('train.pb')
    path('test.pb')
    path('input_genes.txt')
    path('pathway_map.tsv')
    path('model.py')
    path('params.py')
    val(name)

    output:
    path("*history_model.pb"), emit: history
    path("*labels.pb"), emit: labels
    path(name), emit: model
    path('pathways_names.txt'), emit: pathway_labels

    script:
    """
    train_pathway_model.py $name $task.cpus
    """
}
