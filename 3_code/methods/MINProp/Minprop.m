clear;
%refer to : A Heterogeneous Label Progagation Algorithm for Disease Gene Discovery. 
path(path,'../../../2_useful_data');
path(path,'../../common_tool_function');

file_date_time = '2015_8&2016_12';
GP_file_name = ['G_P_network_mappingkey13_' file_date_time '.mat'];

alpha_set = [0.1,0.3,0.5,0.7,0.9];
beta_set = [0.1,0.3,0.5,0.7,0.9];
indistct_parameter_num = 2;

evaluation_indice_set = {'AUC20';'AUC50';'AUC100';'AUC300';'AUC500';'AUC1000';'AUCALL'};
cv_criteria = 'AUC20';
fold_num = 5;
max_ites = 10;
load(GP_file_name,'gene_phenotype_matrix_old', 'gene_phenotype_matrix_newAdded', 'phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','phenotype_id','ppi_matrix');
   
input_parameter_cell = {alpha_set;beta_set;indistct_parameter_num;evaluation_indice_set;cv_criteria;fold_num;max_ites};

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

D = diag(sum(phenotype_similarity_matrix,2));
phenotype_similarity_matrix_norm = D^(-0.5) * phenotype_similarity_matrix * D^(-0.5);

W = Disisolate_ppi(ppi_matrix);
D = diag(sum(W,2));
ppi_matrix_norm = D^(-0.5) * W * D^(-0.5);




matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;...,
    ncbi_gene_id;phenotype_id;matrix_split_output_cell;phenotype_similarity_matrix_norm;ppi_matrix_norm};
matrix_cell_test = {gene_phenotype_matrix_newAdded};
% learning process, using cross-validation to get the best parameter for the
% methods
tic;
[learned_matrix_cell,best_parameter_array,evaluation_parameter_result] = Learn(input_parameter_cell, matrix_cell_train);
toc;

result_file_dir = pwd; %get current directory full name 
%learn_result_cell = {learned_matrix_cell; best_parameter_array; evaluation_parameter_result};
file_key_word = file_date_time;
result_file_name = [result_file_dir '/' 'result' '_' file_key_word '_' datestr(now,30) '.mat' ];  

save(result_file_name,'learned_matrix_cell','best_parameter_array','evaluation_parameter_result');
tic;
    evaluation_index_num = length(evaluation_indice_set);
    evaluation_test = Test(matrix_cell_test,learned_matrix_cell,evaluation_index_num);
toc;
save(result_file_name, 'evaluation_test','-append');