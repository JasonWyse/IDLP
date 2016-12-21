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
    for lambda = alpha_set
        tic;
        max_ites = 10;
        cv_train_parameter_cell = {[lambda]; max_ites;...,
            indistct_num;length(indice_set)};
        [evaluation_result(1:end-indistct_num,i),~] = CrossValidation_Train(cv_train_parameter_cell, matrix_cell_train);
        evaluation_result(end-indistct_num+1:end,i) = [lambda];
        i = i+1;
        toc;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%% cross-validation to choose best parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%
    save('evaluation_result_10.mat','evaluation_result');
    best_parameter_array = Get_best_parameter(evaluation_result,indistct_num,cv_criteria);
    cv_train_parameter_cell = {best_parameter_array; max_ites;indistct_num;length(indice_set)};
    [learned_matrix_cell] = Train2(cv_train_parameter_cell, matrix_cell_train,{gene_phenotype_matrix_old});
   
end

function  [evaluation_result_average,evaluation_result] = CrossValidation_Train(cv_train_parameter_cell, matrix_cell_train)
  %cv_train_parameter_cell = {lambda;eta;gamma; max_ites;...,
   %         indistct_num;length(indice_set)};    
    index_num = cv_train_parameter_cell{end,1};
    fold_num = 5;

    %evaluation_result = zeros(length(initial_matrixFileName_cell),distinct_parameter_num+index_num);
    matrix_split_input_cell = matrix_cell_train(1,1);%if we need multiple matrix to be split, 1-th and 3-th matrix, use cv_train_parameter_cell([1,3],1)
    matrix_split_output_cell = SplitData(matrix_split_input_cell,fold_num);
    %matrix_cell = SplitData(matrix_cell_train,fold_num);
    for i=1:fold_num
        matrix_validation_cell = matrix_split_output_cell(:,i);
        matrix_cell_totrain = MergeData(matrix_split_output_cell,i);   %delete the i-th fold data
        [learned_matrix_cell] = Train2(cv_train_parameter_cell, matrix_cell_train,matrix_cell_totrain);    
        evaluation_result(:,i) = Evaluate(learned_matrix_cell, matrix_validation_cell,index_num);        
    end
    evaluation_result_average = mean(evaluation_result,2);
end