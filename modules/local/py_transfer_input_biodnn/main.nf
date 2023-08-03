// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PY_TRANSFER_INPUT_BIODNN {

    label 'memory_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.4'

    input:
    path('assays.h5')
    path('input_cpgs.pb')
    path('input_cpgs.txt')

    output:
    path("input_list.pb"), emit: input

    script:
    """
    transfer_input_biodnn.py
    """
}
