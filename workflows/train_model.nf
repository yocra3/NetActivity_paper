#!/usr/bin/env nextflow
/*
========================================================================================
Train model on TCGA data
========================================================================================
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////
date = java.time.LocalDate.now()
params.name = ""
params.test_prop = 0.2
params.step = "train"
params.autoencoder = "No"

params.network_params = "a.txt"
params.network = "a.txt"
params.model = "a.txt"
params.config = "a.txt"
params.labels = "a.txt"
params.probes = "a.txt"
params.cpgmap = "a.txt"
params.cpgslist = "a.txt"
params.inputcpgs = "a.txt"
params.gene_mask = "a.txt"
params.checkpoints = "a.txt"

name = params.name
test_prop = params.test_prop
step = params.step
layer = params.layer

ch_input_hdf5 = Channel.fromPath(file("${params.hdf5_file}"))
ch_input_hdf5_se = Channel.fromPath(file("${params.hdf5_se_file}"))

network = file("${params.network}")
network_params = file("${params.network_params}")
model = file("${params.model}")
labels = file("${params.labels}")
probes = file("${params.probes}")
cpgmap = file("${params.cpgmap}")
dnn_structure = file("${params.cpgslist}")
input_cpgs = file("${params.inputcpgs}")
gene_mask = file("${params.gene_mask}")
checkpoints = file("${params.checkpoints}")

random_config = file("${params.config}")
params_name = "v1"
params_name = "${params.params_name}"


include { PY_DIVIDE_TRAIN_TEST } from '../modules/local/py_divide_training_test/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/model_trained/"])
include { RANDOM_SEARCH } from '../modules/local/random_search/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/random_search/"])
include { TRAIN_MODEL } from '../modules/local/py_train_model/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/model_trained/"])
include { PY_EXPORT_RESULTS } from '../modules/local/py_export_results/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/model_trained/"])
include { EXTRACT_MODEL_FEATURES } from '../modules/local/extract_model_features/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/model_features/"])
include { PY_EXTRACT_BIOMODEL_FEATURES } from '../modules/local/py_extract_biomodel_features/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/model_features/"])
include { PY_TRANSFER_INPUT_BIODNN } from '../modules/local/py_transfer_input_biodnn/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}/"])

include { PY_EXPORT_RESULTS_BIO } from '../modules/local/py_export_results_bio/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/model_trained/"])
include { PRUNE_MODEL } from '../modules/local/py_prune_model/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/model_trained/"])


include { TRANSFER_LEARNING } from '../modules/local/transfer_learning/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/transfer_learning/"])
include { PY_EXPORT_RESULTS as PY_EXPORT_RESULTS_TRANSFER} from '../modules/local/py_export_results/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/model_trained/"])
include { PY_EXPORT_RESULTS as PY_EXPORT_RESULTS_AUTOENCOD} from '../modules/local/py_export_results/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/autoencoder_trained/"])
include { PY_EXPORT_RESULTS as PY_EXPORT_RESULTS_AUTOENCOD2} from '../modules/local/py_export_results/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/autoencoder_trained/"])

include { R_FILTER_MISSINGS_HDF5 } from '../modules/local/r_filter_missings_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}"])
include { R_SELECT_AUTOSOMICCPG_HDF5 } from '../modules/local/r_select_autosomic_cpg_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}"])
include { R_FILTER_INVARIANTCPG_HDF5 } from '../modules/local/r_filter_invariantcpg_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}"])
include { R_SUBSTITUTE_MISSING_HDF5 } from '../modules/local/r_substitute_missing_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}"])
include { R_CPGS_MEDIANS } from '../modules/local/r_cpgs_medians/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}"])
include { R_WRITE_PROJECT_LABELS } from '../modules/local/r_write_project_labels/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}"])
include { R_WRITE_INPUT_CPGS } from '../modules/local/r_write_input_cpgs/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}"])
include { R_SELECT_INPUTCPG_HDF5 } from '../modules/local/r_select_inputcpgs_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}"])

include { PY_RESHAPE_HDF5_2D } from '../modules/local/py_reshape_hdf5_2D/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}"])
include { R_CONVERT2D_HDF5 } from '../modules/local/r_convert2D_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}"])
include { PY_RESHAPE_HDF5 } from '../modules/local/py_reshape_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}"])
include { PY_COMBINE_METHY_LABEL_HDF5 } from '../modules/local/py_combine_methy_label_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${params_name}"])
include { TRAIN_AUTOENCODER_MODEL } from '../modules/local/py_train_autoencoder_model/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/autoencoder_trained/"])
include { PY_TRAIN_DNN_BIOLOGICALMODEL } from '../modules/local/py_train_dnn_biologicalmodel/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/model_trained/"])
include { PY_CREATE_GENE_MASK } from '../modules/local/py_create_gene_mask/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/model_trained/"])
include { TRAIN_PATHWAY_MODEL } from '../modules/local/py_train_pathway_model/main.nf' addParams( options: [publish_dir: "${name}/${params_name}/model_trained/"])



workflow  {

  if (step == "preprocessCNN"){
    ch_input = Channel.of("${name}").concat(ch_input_hdf5).concat(ch_input_hdf5_se)
    WF_PREPROCESS_1D(ch_input)
    }
  if (step == "preprocessCNN_trans"){
    ch_input = Channel.of("${name}").concat(ch_input_hdf5).concat(ch_input_hdf5_se)
    WF_PREPROCESS_1D_TRANSF(ch_input, probes)
  }
  if (step == "preprocessDNN"){
    ch_input = Channel.of("${name}").concat(ch_input_hdf5).concat(ch_input_hdf5_se)
    WF_PREPROCESS_DNN(ch_input)
  }
  if (step == "preprocessDNN_trans"){
    ch_input = Channel.of("${name}").concat(ch_input_hdf5).concat(ch_input_hdf5_se)
    WF_PREPROCESS_DNN_TRANSF(ch_input, probes)
  }
  if (step == "preprocess2D"){
    ch_input = Channel.of("${name}").concat(ch_input_hdf5).concat(ch_input_hdf5_se)
    WF_PREPROCESS_2D(ch_input)
  }

  if (step == "random"){
    WF_RANDOM_SEARCH(ch_input_hdf5, test_prop, network, random_config, name)
  }
  if (step == "train"){
    WF_TRAIN_MODEL(ch_input_hdf5, test_prop, network, network_params, name, params.autoencoder)
  }
  if (step == "prune"){
    WF_PRUNE_MODEL(ch_input_hdf5, test_prop, model, network_params, name, params.autoencoder)
  }
  if (step == "train_bio"){
    WF_TRAIN_BIOLOGICALMODEL(ch_input_hdf5, test_prop, network, network_params, probes, cpgmap, name)
  }
  if (step == "train_biopathways"){
    WF_TRAIN_PATHWAY(ch_input_hdf5, test_prop, network, network_params, probes, cpgmap, name, params.autoencoder)
  }
  if (step == "transfer"){
    WF_TRANSFER_MODEL(ch_input_hdf5, test_prop, model, labels, name, layer)
  }
  if (step == "features"){
    EXTRACT_MODEL_FEATURES(ch_input_hdf5, model)
  }
  if (step == "features_bio"){
    ch_input = Channel.of("${name}").concat(ch_input_hdf5).concat(ch_input_hdf5_se)
    WF_EXTRACT_DNN_TRANSF(ch_input,   probes, gene_mask, checkpoints, network, network_params)
  }
  // if (step == "autoencoder"){
  //   //  Old code for autoencoder after CNN
  //   // WF_AUTOENCODER_MODEL(ch_input_hdf5, test_prop, model, labels, name)
  //   WF_AUTOENCODER_MODEL(ch_input_hdf5, test_prop, network, network_params, name)
  // }
}


workflow WF_RANDOM_SEARCH {

  take:
  ch_input_hdf5
  test_prop
  network
  random_config
  name

  main:
  PY_DIVIDE_TRAIN_TEST ( ch_input_hdf5, test_prop )
  RANDOM_SEARCH( PY_DIVIDE_TRAIN_TEST.out.train.collect(), network, random_config, name)

}



workflow WF_TRAIN_BIOLOGICALMODEL {

  take:
  ch_input_hdf5
  test_prop
  network
  network_params
  probes
  cpgmap
  name

  main:
  PY_DIVIDE_TRAIN_TEST ( ch_input_hdf5, test_prop )
  PY_CREATE_GENE_MASK(network_params, probes, cpgmap)
  PY_TRAIN_DNN_BIOLOGICALMODEL( ch_input_hdf5, PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, network, network_params, PY_CREATE_GENE_MASK.out.mask, name )
  PY_EXPORT_RESULTS_BIO( ch_input_hdf5, PY_TRAIN_DNN_BIOLOGICALMODEL.out.history, PY_TRAIN_DNN_BIOLOGICALMODEL.out.labels, PY_DIVIDE_TRAIN_TEST.out.test, PY_CREATE_GENE_MASK.out.mask, PY_TRAIN_DNN_BIOLOGICALMODEL.out.checkpoints, network, network_params, name )

}

workflow WF_TRAIN_PATHWAY {

  take:
  ch_input_hdf5
  test_prop
  network
  network_params
  probes
  pathmap
  name
  autoencoder

  main:
  PY_DIVIDE_TRAIN_TEST ( ch_input_hdf5, test_prop, autoencoder )
  TRAIN_PATHWAY_MODEL( ch_input_hdf5, PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, probes, pathmap, network, network_params, name )
  PY_EXPORT_RESULTS( ch_input_hdf5, TRAIN_PATHWAY_MODEL.out.history, TRAIN_PATHWAY_MODEL.out.model, TRAIN_PATHWAY_MODEL.out.labels, PY_DIVIDE_TRAIN_TEST.out.test, name, autoencoder )

}


workflow WF_PRUNE_MODEL {

  take:
  ch_input_hdf5
  test_prop
  model
  network_params
  name
  autoencoder

  main:
  PY_DIVIDE_TRAIN_TEST ( ch_input_hdf5, test_prop, autoencoder )
  PRUNE_MODEL( ch_input_hdf5, PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, model, network_params, name )
  PY_EXPORT_RESULTS( ch_input_hdf5, PRUNE_MODEL.out.history, PRUNE_MODEL.out.model, PRUNE_MODEL.out.labels, PY_DIVIDE_TRAIN_TEST.out.test, name, autoencoder )

}

workflow WF_TRAIN_MODEL {

  take:
  ch_input_hdf5
  test_prop
  network
  network_params
  name
  autoencoder

  main:
  PY_DIVIDE_TRAIN_TEST ( ch_input_hdf5, test_prop, autoencoder )
  TRAIN_MODEL( ch_input_hdf5, PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, network, network_params, name )
  PY_EXPORT_RESULTS( ch_input_hdf5, TRAIN_MODEL.out.history, TRAIN_MODEL.out.model, TRAIN_MODEL.out.labels, PY_DIVIDE_TRAIN_TEST.out.test, name, autoencoder )

}

workflow WF_TRANSFER_MODEL {

  take:
  ch_input_hdf5
  test_prop
  model
  labels
  name
  layer

  main:
  PY_DIVIDE_TRAIN_TEST ( ch_input_hdf5, test_prop )
  TRANSFER_LEARNING( ch_input_hdf5, PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, model, labels, name, layer )
  PY_EXPORT_RESULTS_TRANSFER( ch_input_hdf5, TRANSFER_LEARNING.out.history, TRANSFER_LEARNING.out.model,  TRANSFER_LEARNING.out.labels, PY_DIVIDE_TRAIN_TEST.out.test, name )

}

workflow WF_AUTOENCODER_MODEL {

  take:
  ch_input_hdf5
  test_prop
  model
  params
  // labels
  name
  //  Old code for autoencoder after CNN

  main:
  PY_DIVIDE_TRAIN_TEST ( ch_input_hdf5, test_prop )
  TRAIN_AUTOENCODER_MODEL( PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, model, params, name )

  // TRAIN_AUTOENCODER_MODEL( PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, model, name )
  // PY_EXPORT_RESULTS_AUTOENCOD( ch_input_hdf5, TRAIN_AUTOENCODER_MODEL.out.old_history, TRAIN_AUTOENCODER_MODEL.out.model_old,  labels, PY_DIVIDE_TRAIN_TEST.out.test, name )

}


workflow WF_PREPROCESS_1D {

  take:
  ch_input_hdf5


  main:
  R_FILTER_MISSINGS_HDF5( ch_input_hdf5.collect() )
  R_SUBSTITUTE_MISSING_HDF5( R_FILTER_MISSINGS_HDF5.out.res )

  R_CPGS_MEDIANS( R_FILTER_MISSINGS_HDF5.out.res )
  R_WRITE_PROJECT_LABELS( ch_input_hdf5.collect() )

  PY_RESHAPE_HDF5( R_SUBSTITUTE_MISSING_HDF5.out.h5, R_WRITE_PROJECT_LABELS.out )

}

workflow WF_PREPROCESS_DNN {

  take:
  ch_input_hdf5


  main:
  R_FILTER_MISSINGS_HDF5( ch_input_hdf5.collect() )
  R_SELECT_AUTOSOMICCPG_HDF5( R_FILTER_MISSINGS_HDF5.out.res )
  R_SUBSTITUTE_MISSING_HDF5( R_SELECT_AUTOSOMICCPG_HDF5.out.res )

  R_CPGS_MEDIANS( R_SELECT_AUTOSOMICCPG_HDF5.out.res )
  R_WRITE_INPUT_CPGS( R_SELECT_AUTOSOMICCPG_HDF5.out.res )
  R_WRITE_PROJECT_LABELS( ch_input_hdf5.collect() )

  PY_COMBINE_METHY_LABEL_HDF5( R_SUBSTITUTE_MISSING_HDF5.out.h5, R_WRITE_PROJECT_LABELS.out )

}

workflow WF_PREPROCESS_2D {

  take:
  ch_input_hdf5


  main:
  R_FILTER_MISSINGS_HDF5( ch_input_hdf5.collect() )
  R_SUBSTITUTE_MISSING_HDF5( R_FILTER_MISSINGS_HDF5.out.res )

  R_CPGS_MEDIANS( R_FILTER_MISSINGS_HDF5.out.res )
  R_WRITE_PROJECT_LABELS( ch_input_hdf5.collect() )
  R_CONVERT2D_HDF5(R_SUBSTITUTE_MISSING_HDF5.out.res )
  PY_RESHAPE_HDF5_2D( R_CONVERT2D_HDF5.out.res, R_WRITE_PROJECT_LABELS.out )


}

workflow WF_PREPROCESS_1D_TRANSF {

  take:
  ch_input_hdf5
  probes

  main:
  R_SUBSTITUTE_MISSING_HDF5( ch_input_hdf5.collect() )
  R_SELECT_INPUTCPG_HDF5( R_SUBSTITUTE_MISSING_HDF5.out.res , probes )

  R_WRITE_PROJECT_LABELS( ch_input_hdf5.collect() )

  PY_RESHAPE_HDF5( R_SELECT_INPUTCPG_HDF5.out.h5, R_WRITE_PROJECT_LABELS.out )

}

workflow WF_PREPROCESS_DNN_TRANSF {

  take:
  ch_input_hdf5
  probes

  main:
  R_SUBSTITUTE_MISSING_HDF5( ch_input_hdf5.collect() )
  R_SELECT_INPUTCPG_HDF5( R_SUBSTITUTE_MISSING_HDF5.out.res , probes )

  R_WRITE_PROJECT_LABELS( ch_input_hdf5.collect() )

  PY_COMBINE_METHY_LABEL_HDF5( R_SELECT_INPUTCPG_HDF5.out.h5, R_WRITE_PROJECT_LABELS.out )

}


workflow WF_EXTRACT_DNN_TRANSF {

  take:
  ch_input_hdf5
  cpg_median
  gene_mask
  checkpoints
  network
  network_params

  main:

  R_SUBSTITUTE_MISSING_HDF5( ch_input_hdf5.collect() )
  R_SELECT_INPUTCPG_HDF5( R_SUBSTITUTE_MISSING_HDF5.out.res , cpg_median )
  PY_EXTRACT_BIOMODEL_FEATURES(R_SELECT_INPUTCPG_HDF5.out.h5, gene_mask, checkpoints, network, network_params)

}
