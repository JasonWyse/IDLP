<<<<<<< HEAD
function [learned_matrix_cell, best_parameter_array, evaluation_result] = Learn(input_parameter_cell, matrix_cell_train, ...,
    initial_matrixFileName_cell)
    %%%%%%%%%%%%%%%%%%%%%%%%%% parsing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha_set = input_parameter_cell{1,1};
    indice_set = input_parameter_cell{2,1};
    cv_criteria =  input_parameter_cell{3,1};
    gene_phenotype_matrix_old = matrix_cell_train{1,1};
    %%%%%%%%%%%%%%%%%%%%%%%%%% parsing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    indistct_num = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%% cross-validation to choose best parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%
    evaluation_result = zeros(indistct_num+length(indice_set),length(alpha_set));
    i = 1;
    for alpha = alpha_set
        tic;
        max_ites = 20;
        cv_train_parameter_cell = {alpha; max_ites;indistct_num;length(indice_set)};
        [evaluation_result(1:end-indistct_num,i),~] = CrossValidation_Train(cv_train_parameter_cell, matrix_cell_train);
        evaluation_result(end-indistct_num+1:end,i) = alpha;
=======
function [learned_matrix_cell, best_parameter_array, evaluation_result] = Learn(input_parameter_cell, matrix_cell_train)
    %%%%%%%%%%%%%%%%%%%%%%%%%% parsing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha_set = input_parameter_cell{1,1};
    indistct_parameter_num = input_parameter_cell{2,1};
    evaluation_indice_set = input_parameter_cell{3,1};
    cv_criteria =  input_parameter_cell{4,1};
    fold_num = input_parameter_cell{5,1};
    max_ites = input_parameter_cell{6,1};
    gene_phenotype_matrix_old = matrix_cell_train{1,1};
    %%%%%%%%%%%%%%%%%%%%%%%%%% parsing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% cross-validation to choose best parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%   
    evaluation_index_num = length(evaluation_indice_set);
    evaluation_result = zeros(indistct_parameter_num+evaluation_index_num,length(alpha_set));
    i = 1;    
    for alpha = alpha_set
        tic;        
        cv_train_parameter_cell = {alpha; max_ites;indistct_parameter_num;evaluation_index_num;fold_num};
        [evaluation_result(1:end-indistct_parameter_num,i),~] = CrossValidation_Train(cv_train_parameter_cell, matrix_cell_train);
        evaluation_result(end-indistct_parameter_num+1:end,i) = alpha;
>>>>>>> master
        i = i+1;
        toc;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%% cross-validation to choose best parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%
    save('evaluation_result.mat','evaluation_result');
<<<<<<< HEAD
    best_parameter_array = Get_best_parameter(evaluation_result,indistct_num,cv_criteria);
    cv_train_parameter_cell = {best_parameter_array; max_ites;indistct_num;length(indice_set)};
    [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train,gene_phenotype_matrix_old);
   
end

function  [evaluation_result_average,evaluation_result] = CrossValidation_Train(cv_train_parameter_cell, matrix_cell_train)
    index_num = cv_train_parameter_cell{end,1};
    fold_num = 2;
    alpha = cv_train_parameter_cell{1,1};
    %evaluation_result = zeros(length(initial_matrixFileName_cell),distinct_parameter_num+index_num);
    matrix_split_input_cell = matrix_cell_train(1,1);%if we need multiple matrix to be split, 1-th and 3-th matrix, use cv_train_parameter_cell([1,3],1)
    matrix_split_output_cell = SplitData(matrix_split_input_cell,fold_num);
    %matrix_cell = SplitData(matrix_cell_train,fold_num);
    for i=1:fold_num
        matrix_validation_cell = matrix_split_output_cell(:,1);
        matrix_cell_totrain = MergeData(matrix_split_output_cell,i);   %delete the i-th fold data
        [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train,matrix_cell_totrain);    
        evaluation_result(:,i) = Evaluate(learned_matrix_cell, matrix_validation_cell,index_num);        
    end
    evaluation_result_average = mean(evaluation_result,2);
=======
    best_parameter_array = Get_best_parameter(evaluation_result,indistct_parameter_num,cv_criteria);
    cv_train_parameter_cell = {best_parameter_array; max_ites;indistct_parameter_num;evaluation_index_num};
    [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train,{gene_phenotype_matrix_old},fold_num+1);
   
end

function  [evaluation_result_average,evaluation_folds_result] = CrossValidation_Train(cv_train_parameter_cell, matrix_cell_train)
    evaluation_index_num = cv_train_parameter_cell{end-1,1};
    fold_num = cv_train_parameter_cell{end,1};    
    %evaluation_result = zeros(length(initial_matrixFileName_cell),distinct_parameter_num+index_num);
    %matrix_split_input_cell = matrix_cell_train(1,1);%if we need multiple matrix to be split, 1-th and 3-th matrix, use cv_train_parameter_cell([1,3],1)
    %matrix_split_output_cell = SplitData(matrix_split_input_cell,fold_num);
    matrix_split_output_cell = matrix_cell_train{6,1};
    %matrix_cell = SplitData(matrix_cell_train,fold_num);
    evaluation_folds_result = zeros(evaluation_index_num,fold_num);
    for i=1:fold_num
        matrix_oneFold_cell_validate = matrix_split_output_cell(:,i);
        matrix_mergeFolds_cell_train = MergeData(matrix_split_output_cell,i);   %delete the i-th fold data
        [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train,matrix_mergeFolds_cell_train,i);    
        evaluation_folds_result(:,i) = Evaluate(learned_matrix_cell, matrix_oneFold_cell_validate,evaluation_index_num);        
    end
    evaluation_result_average = mean(evaluation_folds_result,2);
>>>>>>> master
end