clear;
%refer to : Walking the Interactome for Prioritization of Candidate Disease Genes
path(path,'../../../2_useful_data');
path(path,'../../common_tool_function');
file_date_time = '2015_8&2016_12';

ppi_key_word = 'exist_old_';
GP_file_name = ['G_P_network_' ppi_key_word 'mappingkey13_' file_date_time '.mat'];

beta_set = [0.1,0.3,0.5,0.7,0.9];
indice_set = {'AUC20';'AUC50';'AUC100';'AUC300';'AUC500';'AUC1000';'AUCALL'};
indice_num = 7;
indistct_parameter_num = 1;
fold_num = 5;
cv_criteria = 'AUC20';
load(GP_file_name,'gene_phenotype_matrix_old', 'gene_phenotype_matrix_newAdded', 'phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','phenotype_id','ppi_matrix');
   
input_parameter_cell = {beta_set;indistct_parameter_num;indice_set;cv_criteria;fold_num};
matrix_split_input_cell = {gene_phenotype_matrix_old};
%-->gene-phenotype associations are divided into X folds, these folds are the same for different  
%-->parameters, we can save these folds into a mat file in order to avoid
%-->duplication of gene-phenotype association division. It can greatly
%-->shorten the execute time
useful_data_dir = '../../../2_useful_data/';
filename = [useful_data_dir 'matrix_split_output_cell_' num2str(fold_num) 'folds.mat'];
if exist(filename,'file')
    load(filename,'matrix_split_output_cell');
else    
    [ matrix_split_output_cell ] = SplitData( matrix_split_input_cell,fold_num );
    %Y_hat = Initialize_Y_prince(gene_phenotype_matrix, phenotype_similarity_matrix);    
    save(filename,'matrix_split_output_cell');
end

array = sum(ppi_matrix,2);
ppi_matrix_diag = diag(array);
ppi_matrix_minus = ppi_matrix_diag - ppi_matrix;

matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;...,
    ppi_matrix;ncbi_gene_id;phenotype_id;matrix_split_output_cell;ppi_matrix_minus};
matrix_cell_test = {gene_phenotype_matrix_newAdded};
% learning process, using cross-validation to get the best parameter for the
% methods
tic;
[learned_matrix_cell,best_parameter_array,evaluation_result] = Learn(input_parameter_cell, matrix_cell_train);
toc;
result_file_dir = pwd; %get current directory full name 
file_key_word = file_date_time;
result_file_name = [result_file_dir '/' 'result' '_' ppi_key_word file_key_word '_' datestr(now,30) '.mat' ];  

save(result_file_name,'learned_matrix_cell','best_parameter_array','evaluation_result');
tic;
    evaluation_test = Test(matrix_cell_test,learned_matrix_cell,indice_num);
toc;
save(result_file_name, 'evaluation_test','-append');