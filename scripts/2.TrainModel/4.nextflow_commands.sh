#'#################################################################################
#'#################################################################################
# Nextflow commands for applying NetActivity
#'#################################################################################
#'#################################################################################

# Train main model in GTEx
## Train with all GO terms -> full training
## Train GTEx gexp protein coding pathways  standardized - v3.11
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_full_v3.11 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -resume -profile docker


## GTEx gexp pathways - v3.11 features
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_full_v3.11 --step features \
--model  results/GTEx_coding/paths_all_full_v3.11/model_trained/GTEx_coding -resume -profile docker

## Train TCGA gexp pathways - v3.11 - 5 initializations
for i in {a..e}
do
echo paths_all_full_v3.11${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_full_v3.11${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_full_v3.11${i}  --step features \
--model  results/GTEx_coding/paths_all_full_v3.11${i}/model_trained/GTEx_coding -profile docker
done

## Train with subset of pathways 1
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt1_full_v3.11 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt_gene_map.tsv -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt1_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt1_full_v3.11/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo paths_filt1_full_v3.11${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt1_full_v3.11${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt1_full_v3.11${i}  --step features \
--model  results/GTEx_coding/paths_filt1_full_v3.11${i}/model_trained/GTEx_coding -profile docker
done

## Train with subset of pathways 2
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_v3.11 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo paths_filt2_full_v3.11${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_v3.11${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_v3.11${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11${i}/model_trained/GTEx_coding -profile docker
done


# Train model in GTEx with different configurations
## Train GTEx all paths only step 1
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.8 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.8 --step features \
--model  results/GTEx_coding/paths_all_pretrain_v3.8/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.8${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.8${i}  --step features \
--model  results/GTEx_coding/paths_all_pretrain_v3.8${i}/model_trained/GTEx_coding -profile docker
done


# Train GTEx selected paths only step 1
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_pre_v3.8 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_pre_v3.8 --step features \
--model  results/GTEx_coding/paths_filt2_pre_v3.8/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo paths_filt2_pre_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_pre_v3.8${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_pre_v3.8${i}  --step features \
--model  results/GTEx_coding/paths_filt2_pre_v3.8${i}/model_trained/GTEx_coding -profile docker
done

# Train GTEx all paths step 1 + step 3
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.10 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v10.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.10 --step features \
--model  results/GTEx_coding/paths_all_pretrain_v3.10/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v3.10${i}
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.10${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v10.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.10${i}  --step features \
--model  results/GTEx_coding/paths_all_pretrain_v3.10${i}/model_trained/GTEx_coding -profile docker
done


# Train GTEx selected paths step 1 + step 3
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_unfrozen_v3.10 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v10.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_unfrozen_v3.10 --step features \
--model  results/GTEx_coding/paths_filt2_unfrozen_v3.10/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo unfrozen ${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_unfrozen_v3.10${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v10.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_unfrozen_v3.10${i}  --step features \
--model  results/GTEx_coding/paths_filt2_unfrozen_v3.10${i}/model_trained/GTEx_coding -profile docker
done


# Train GTEx selected paths only step 3
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_noprime_v3.12 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v12.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_noprime_v3.12 --step features \
--model  results/GTEx_coding/paths_filt2_full_noprime_v3.12/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_noprime_v3.12${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v12.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_noprime_v3.12${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_noprime_v3.12${i}/model_trained/GTEx_coding -profile docker
done


# Train GTEx selected paths Dropout - no step 2
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_noprime_v3.7 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_noprime_v3.7 --step features \
--model  results/GTEx_coding/paths_filt2_full_drop_noprime_v3.7/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_noprime_v3.7${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_noprime_v3.7${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_drop_noprime_v3.7${i}/model_trained/GTEx_coding -profile docker
done

# Train GTEx selected paths dropout + step 2
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_prime_v3.9 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v9.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_prime_v3.9 --step features \
--model  results/GTEx_coding/paths_filt2_full_drop_prime_v3.9/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_prime_v3.9${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v9.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_prime_v3.9${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_drop_prime_v3.9${i}/model_trained/GTEx_coding -profile docker
done


# Train GTEx selected paths - 2-layer decoder
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_postdense_v4.3 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_postdense_v4.3 --step features \
--model  results/GTEx_coding/paths_filt2_full_postdense_v4.3/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v4.3${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_postdense_v4.3${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_postdense_v4.3${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_postdense_v4.3${i}/model_trained/GTEx_coding -profile docker
done


# Train GTEx selected paths - 2-layer encoder and 2-layer decoder
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_prepostdense_v5.3 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_prepostdense_v5.3 --step features \
--model  results/GTEx_coding/paths_filt2_full_prepostdense_v5.3/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo prepostdense_v5.3${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_prepostdense_v5.3${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker


nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_prepostdense_v5.3${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_prepostdense_v5.3${i}/model_trained/GTEx_coding -profile docker
done


# Train GTEx selected paths - 2-layer encoder
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_predense_v6.2 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network6_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network6.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_predense_v6.2 --step features \
--model  results/GTEx_coding/paths_filt2_full_predense_v6.2/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_predense_v6.2${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network6_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network6.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv -resume -profile docker


nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_predense_v6.2${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_predense_v6.2${i}/model_trained/GTEx_coding -profile docker
done


# Train subsets of gene sets
sizes=(30 100 300 1000)

for size in "${sizes[@]}"
do
  echo "Size: $size"
  for rep in {1..10}
  do
    for i in {a..f}
    do
        nextflow run yocra3/NetActivityTrain --data_prefix results/GTEx/vst_all_group_ --gene_mask results/subsetGeneSets/subset_${rep}_N${size}_geneMap.tsv \
        --network conf/network_params/dnn_gexp_autoencod_pathway_network3.py \
        --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py --outdir results/subsetGeneSets/models/size_${size}/${rep}/${i}/ -profile docker
    done
  done
done

sizes=(30 100 300 1000)

for size in "${sizes[@]}"
do
  echo "Size: $size"
  for rep in {1..10}
  do
         nextflow run yocra3/NetActivityTrain --data_prefix results/GTEx/vst_all_group_ --gene_mask results/subsetGeneSets/subset_${rep}_N${size}_geneMap.tsv \
        --network conf/network_params/dnn_gexp_autoencod_pathway_network3.py \
        --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py --outdir results/subsetGeneSets/models/size_${size}/${rep}/f/ -profile docker
  done
done


# Train random gene sets
for i in {a..f}
do
  nextflow run yocra3/NetActivityTrain --data_prefix results/GTEx/vst_all_group_ --gene_mask results/randomGeneSets/random_geneMap.tsv \
  --network conf/network_params/dnn_gexp_autoencod_pathway_network3.py \
  --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py --outdir results/randomGeneSets/${i}/ -profile docker
done


# Train model in TCGA data
## TCGA gexp all - v3.11 features
nextflow run yocra3/NetActivityTrain --data_prefix results/TCGA_gexp_combat_coding/vsd_norm_group --gene_mask results/TCGA_gexp_combat_coding/gene_map.tsv \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py \
--network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py --outdir results/TCGA_model/ -profile docker


#


# Compute GSAS in different datasets
## Compute GSAS in PRAD
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/prad_assay_reshaped_standardized.h5  \
--name GTEx_coding_PRAD --params_name paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker

## Compute GSAS in GSE169038
nextflow run workflows/train_model.nf --hdf5_file results/GSE169038/prad_array_reshaped_standardized.h5  \
--name GSE169038 --params_name paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker

## Compute GSAS in SRP042228
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_coding_std_gse.h5  \
--name SRP042228 --params_name paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker


