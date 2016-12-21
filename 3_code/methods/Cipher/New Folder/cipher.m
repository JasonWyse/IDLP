clear;

path(path,'../../../2_useful_data');
path(path,'../../common_tool_function');
GP_file_name = 'G_P_network_mappingkey13_2015_8&2016_12.mat';

method_str = 'sp';  %sp  dn
indice_set = {'AUC50';'AUC100';'AUC300';'AUC500';'AUC1000';'AUCALL'};
cv_criteria = 'AUC50';
load(GP_file_name,'gene_phenotype_matrix_old', 'gene_phenotype_matrix_newAdded', 'phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','phenotype_id','ppi_matrix');
ppi_sp = graphallshortestpaths(sparse(ppi_matrix));

input_parameter_cell = {method_str;indice_set;cv_criteria};
matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ppi_sp;ncbi_gene_id;phenotype_id};
matrix_cell_test = {gene_phenotype_matrix_newAdded};
% learning process, using cross-validation to get the best parameter for the
% methods
tic;
[learned_matrix_cell,best_parameter_array,evaluation_result] = Learn(input_parameter_cell, matrix_cell_train);
toc;
save('result.mat','learned_matrix_cell','best_parameter_array');
tic;
    evaluation_test = Test(matrix_cell_test,learned_matrix_cell,6);
toc;
save('result.mat', 'evaluation_test','-append');