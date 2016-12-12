clear;
%refer to : Genome-wide inferring gene¨Cphenotype relationship by walking on the heterogeneous network.pdf
path(path,'../../../2_useful_data');
path(path,'../../common_tool_function');
GP_file_name = 'G_P_network_mappingkey13_2015_8&2016_12.mat';

lambda_set = [0.5];
eta_set = [0.3,0.5];
gamma_set = [0.7];
indice_set = {'AUC50';'AUC100';'AUC300';'AUC500';'AUC1000';'AUCALL'};
cv_criteria = 'AUC50';
load(GP_file_name,'gene_phenotype_matrix_old', 'gene_phenotype_matrix_newAdded', 'phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','phenotype_id','ppi_matrix');
   
input_parameter_cell = {lambda_set;eta_set;gamma_set;indice_set;cv_criteria};
matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ncbi_gene_id;phenotype_id};
matrix_cell_test = {gene_phenotype_matrix_newAdded};
% learning process, using cross-validation to get the best parameter for the
% methods
tic;
[learned_matrix_cell,best_parameter_array,evaluation_result] = Learn(input_parameter_cell, matrix_cell_train);
toc;
save('train.mat','learned_matrix_cell','best_parameter_array','evaluation_result');
tic;
    evaluation_test = Test(matrix_cell_test,learned_matrix_cell);
toc;
 