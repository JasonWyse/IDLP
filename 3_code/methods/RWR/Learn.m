function [learned_matrix_cell, best_parameter_array, evaluation_result] = Learn(input_parameter_cell, matrix_cell_train)
    %%%%%%%%%%%%%%%%%%%%%%%%%% parsing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda_set = input_parameter_cell{1,1};
    indistct_parameter_num = input_parameter_cell{2,1};
    evaluation_indice_set = input_parameter_cell{3,1};
    cv_criteria =  input_parameter_cell{4,1};
    fold_num = input_parameter_cell{5,1};
    max_ites = input_parameter_cell{6,1};
    gene_phenotype_matrix_old = matrix_cell_train{1,1};
    %%%%%%%%%%%%%%%%%%%%%%%%%% parsing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% cross-validation to choose best parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%
    evaluation_index_num = length(evaluation_indice_set);
    evaluation_result = zeros(indistct_parameter_num+evaluation_index_num,length(lambda_set));
    i = 1;
    for lambda = lambda_set
        tic;
        cv_train_parameter_cell = {lambda; max_ites;indistct_parameter_num;evaluation_index_num;fold_num};
        [evaluation_result(1:end-indistct_parameter_num,i),~] = CrossValidation_Train(cv_train_parameter_cell, matrix_cell_train);
        evaluation_result(end-indistct_parameter_num+1:end,i) = lambda;
        i = i+1;
        toc;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%% cross-validation to choose best parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%
    best_parameter_array = Get_best_parameter(evaluation_result,indistct_parameter_num,cv_criteria);
    cv_train_parameter_cell = {best_parameter_array; max_ites;indistct_parameter_num;evaluation_index_num};
    [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train,fold_num+1);
   
end

function  [evaluation_result_average,evaluation_result] = CrossValidation_Train(cv_train_parameter_cell, matrix_cell_train)
    evaluation_index_num = cv_train_parameter_cell{end-1,1};
    fold_num = cv_train_parameter_cell{end,1};   
    matrix_split_output_cell = matrix_cell_train{6,1};
    evaluation_result = zeros(evaluation_index_num,fold_num);
    for i=1:fold_num
        matrix_validation_cell = matrix_split_output_cell(:,i);
        [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train,i);    
        evaluation_result(:,i) = Evaluate(learned_matrix_cell, matrix_validation_cell,evaluation_index_num);        
    end
    evaluation_result_average = mean(evaluation_result,2);
end