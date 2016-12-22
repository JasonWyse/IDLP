function [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train,fold_idx)
    %cv_train_parameter_cell = {alpha; max_ites;indistct_num;length(indice_set)};
    %matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ncbi_gene_id;phenotype_id};
    %matrix_cell_totrain = MergeData(matrix_split_output_cell,i);
    parameter = cv_train_parameter_cell{1,1};
    lambda = parameter(1,1);
    max_ites = cv_train_parameter_cell{2,1};
    
    
    S =  matrix_cell_train{7,1};
    Y_hat_folds_cell = matrix_cell_train{8,1};

%     array = 1./sum(gene_phenotype_matrix_old);
%     array(isinf(array)) = 0;
%     Y = gene_phenotype_matrix_old* diag(array);
    Y = Y_hat_folds_cell{1,fold_idx};
    F_before = Y;
    for step = 1:max_ites
      F_now = (1-lambda) * S * F_before + lambda * Y; 
      F_before = F_now;
    end
    gene_phenotype_score_matrix = F_now;


    learned_matrix_cell = {gene_phenotype_score_matrix};
end

