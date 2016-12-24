clear;

path(path,'../../../2_useful_data');
file_date_time = '2015_8&2016_12';
key_word = 'exist_old_';
GP_file_name = ['G_P_network_' key_word 'mappingkey13_' file_date_time '.mat'];

method_str = 'sp';  %sp  dn
evaluation_indice_set = {'AUC20';'AUC50';'AUC100';'AUC300';'AUC500';'AUC1000';'AUCALL'};
load(GP_file_name,'gene_phenotype_matrix_old', 'gene_phenotype_matrix_newAdded', 'phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','phenotype_id','ppi_matrix','gene_phenotype_matrix_new');
ppi_sp = graphallshortestpaths(sparse(ppi_matrix));
if strcmp(method_str,'sp') == 1
    dm = 1;
    fn = 'cipher_result_test_sp.mat';
elseif strcmp(method_str,'dn') == 1
    dm =2;
    fn = 'cipher_result_test_dn.mat';
end

[~,cols] = size(gene_phenotype_matrix_old);

[allGene_num, ~] = size(gene_phenotype_matrix_old);
allgene_profile = zeros(size(gene_phenotype_matrix_old'));
for i = 1:allGene_num 
    allgene_profile(:,i) = get_gene_profile(i, gene_phenotype_matrix_old, ppi_sp, dm);
end 
matrix_conbine = [allgene_profile,phenotype_similarity_matrix(:,1:cols)];
pearson_matrix = corrcoef(matrix_conbine);
gene_phenotype_score_matrix = pearson_matrix(1:allGene_num,allGene_num+1:end); 


cROC = zeros(1,7);
[~,avgROC] = ROC_main(gene_phenotype_matrix_newAdded',gene_phenotype_score_matrix');
cROC(1,:) = avgROC;

result_file_dir = pwd;
result_file_name = [result_file_dir '/' 'result' '_' key_word file_date_time '_' datestr(now,30) '_' method_str '_' '.mat' ];  


save (result_file_name,'gene_phenotype_score_matrix','cROC');