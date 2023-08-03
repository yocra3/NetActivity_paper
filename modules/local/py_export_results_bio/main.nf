// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PY_EXPORT_RESULTS_BIO {

    label 'memory_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.3b'

    input:
    path('assay_reshaped.h5')
    path('history_model.pb')
    path('labels.pb')
    path('test.pb')
    path('gene_mask.pb')
    path(checkpoints)
    path('model.py')
    path('params.py')
    val(name)

    output:
    path("*.tsv"), emit: tsv
    path("*.txt"), emit: txt

    script:
    """
    export_results_bio.py $name
    """
}
