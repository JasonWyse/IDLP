clear;
%refer to : Walking the Interactome for Prioritization of Candidate Disease Genes
path(path,'../../../2_useful_data');
path(path,'../../common_tool_function');
GP_file_name = 'G_P_network_mappingkey13_2015_8&2016_12.mat';

beta_set = [0.1,0.3,0.5,0.7,0.9];
indice_set = {'AUC50';'AUC100';'AUC300';'AUC500';'AUC1000';'AUCALL'};
indice_num = 6;
cv_criteria = 'AUC50';
load(GP_file_name,'gene_phenotype_matrix_old', 'gene_phenotype_matrix_newAdded', 'phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','phenotype_id','ppi_matrix');
   
input_parameter_cell = {beta_set;indice_set;cv_criteria};
matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ncbi_gene_id;phenotype_id};
matrix_cell_test = {gene_phenotype_matrix_newAdded};
% learning process, using cross-validation to get the best parameter for the
% methods
tic;ye
[learned_matrix_cell,best_parameter_array,evaluation_result] = Learn(input_parameter_cell, matrix_cell_train);
toc;
save('result.mat','learned_matrix_cell','best_parameter_array');
tic;
    evaluation_test = Test(matrix_cell_test,learned_matrix_cell,indice_num);
toc;
save('result.mat', 'evaluation_test','-append');