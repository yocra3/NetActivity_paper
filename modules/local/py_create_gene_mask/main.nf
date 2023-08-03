// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PY_CREATE_GENE_MASK {

    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.3b'

    input:
    path('params.py')
    path('inputcpgs.txt')
    path('cpgs_map.txt')

    output:
    path("gene_mask.pb"), emit: mask

    script:
    """
    create_gene_mask.py
    """
}
