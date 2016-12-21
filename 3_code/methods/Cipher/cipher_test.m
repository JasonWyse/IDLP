function cipher_test(dist_metric)

if strcmp(dist_metric,'sp') == 1
    dm = 1;
    fn = 'cipher_result_test_sp.mat';
elseif strcmp(dist_metric,'dn') == 1
    dm = 2;
    fn = 'cipher_result_test_dn.mat';
end

GP_file_name = 'G_P_network_mappingkey13_2015_8&2016_12.mat';
load(GP_file_name,'gene_phenotype_matrix_old', 'gene_phenotype_matrix_newAdded', 'phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','phenotype_id','ppi_matrix');
   
[rows, cols] = size(gene_phenotype_matrix_old);
ppi_sp = graphallshortestpaths(sparse(ppi_matrix));
gene_phenotype_score_matrix = zeros(rows, cols);
gene_phenotype_score_matrix = CipherRank(phenotype_similarity_matrix, gene_phenotype_matrix_old, 1, cols, ppi_sp, dm);

cROC = zeros(1,6);
[~,avgROC] = ROC_main(gene_phenotype_matrix_newAdded',gene_phenotype_score_matrix');
cROC(1,:) = avgROC;

save (fn,'gene_phenotype_score_matrix','cROC');
end
