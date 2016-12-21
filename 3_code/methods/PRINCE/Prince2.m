clear;
%refer to : Associating genes and Protein Complexes with Disease via Network Progagation
path(path,'../../../2_useful_data');
path(path,'../../common_tool_function');

%%%%%%%%%%%%%%%%%%%% load data file %%%%%%%%%%%%%%%%%%%%%%
file_date_time = '2015_8&2016_12';
GP_file_name = ['G_P_network_mappingkey13_' file_date_time '.mat'];
load(GP_file_name,'gene_phenotype_matrix_old', 'gene_phenotype_matrix_newAdded', 'phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','phenotype_id','ppi_matrix');
%%%%%%%%%%%%%%%%%%%% load data file %%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%% input parameter assignment %%%%%%%%%%%%%%%
alpha_set = [0.1];
indistct_parameter_num = 1;
evaluation_indice_set = {'AUC20';'AUC50';'AUC100';'AUC300';'AUC500';'AUC1000';'AUCALL'};
%index_cell = {'AUC50';'AUC100';'AUC300'; 'AUC500';'AUC1000';'AUCALL'};
cv_criteria = 'AUC20';
fold_num = 5;
max_ites = 20;
%%%%%%%%%%%%%%%%%%%% input parameter assignment %%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% prepare data cell %%%%%%%%%%%%%%%%%%%
input_parameter_cell = {alpha_set;indistct_parameter_num;evaluation_indice_set;cv_criteria;fold_num;max_ites};
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

filename = 'Y_hat_folds_cell.mat';% seed disease genes don't change for each query disease, 
                                  % it can be prepared in advance to
                                  % shorten the execute time
if exist(filename,'file')
    load('Y_hat_folds_cell.mat','Y_hat_folds_cell');
else    
    [ Y_hat_folds_cell ] = Get_Y_hat_folds( matrix_split_output_cell, phenotype_similarity_matrix );
    Y_hat_folds_cell{1,fold_num+1} = Get_Y_hat_allFolds( gene_phenotype_matrix_old, phenotype_similarity_matrix );    
    save(filename,'Y_hat_folds_cell');
end
matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;...,
    ncbi_gene_id;phenotype_id;matrix_split_output_cell;Y_hat_folds_cell};
matrix_cell_test = {gene_phenotype_matrix_newAdded};
%%%%%%%%%%%%%%%%%%%% prepare data cell %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% learning process %%%%%%%%%%%%%%%%%%%%%%
% learning process, using cross-validation to get the best parameter for the methods
tic;
[learned_matrix_cell,best_parameter_array,evaluation_parameter_result] = Learn(input_parameter_cell, matrix_cell_train);
toc;
result_file_dir = pwd; %get current directory full name 
%learn_result_cell = {learned_matrix_cell; best_parameter_array; evaluation_parameter_result};
file_key_word = file_date_time;
result_file_name = [result_file_dir '/' 'result' '_' file_key_word '_' datestr(now,30) '.mat' ];  
save(result_file_name, 'learned_matrix_cell', 'best_parameter_array', 'evaluation_parameter_result',...,
    'max_ites');
%save('result_logic.mat','learned_matrix_cell','best_parameter_array');
%%%%%%%%%%%%%%%%%%%% learning process %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% test process %%%%%%%%%%%%%%%%%%%%%%
tic;
    evaluation_index_num = length(evaluation_indice_set);
    evaluation_test = Test(matrix_cell_test,learned_matrix_cell,evaluation_index_num);
toc;
save(result_file_name, 'evaluation_test','-append');
%save('result_logic.mat', 'evaluation_test','-append');
%%%%%%%%%%%%%%%%%%%% test process %%%%%%%%%%%%%%%%%%%%%%