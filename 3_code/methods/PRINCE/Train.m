function [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train,...,
    matrix_cell_totrain,fold_idx)
    %cv_train_parameter_cell = {alpha; max_ites;indistct_num;length(indice_set)};
    %matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ncbi_gene_id;phenotype_id};
    %matrix_cell_totrain = MergeData(matrix_split_output_cell,i);
    parameter = cv_train_parameter_cell{1,1};
    alpha = parameter(1,1);
    max_ites = cv_train_parameter_cell{2,1};    
    %gene_phenotype_matrix_old = matrix_cell_totrain{1,1};
%    matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;...,
%    ncbi_gene_id;phenotype_id;matrix_split_output_cell};
%    phenotype_similarity_matrix = matrix_cell_train{2,1};
    ppi_matrix = matrix_cell_train{3,1};  
    Y_hat_folds_cell = matrix_cell_train{7,1};
    
    W = Disisolate_ppi(ppi_matrix);
    D = diag(1./sqrt(sum(W,2)));
    S = D*W*D;   
%     D = diag(sum(W,2));
%     S = D^(-0.5) * W * D^(-0.5);
    % phenotype_similarity_matrix = 1./(1+exp((-15*phenotype_similarity_matrix)+log(9999)));
     Y_temp = Y_hat_folds_cell{1,fold_idx};
     Y = 1./(1+exp((-15*Y_temp)+log(9999)));
    
    F_before = Y;
    for step = 1:max_ites
      F_now = alpha * S * F_before + (1-alpha) * Y; 
      F_before = F_now;
    end
    gene_phenotype_score_matrix = F_now;
    learned_matrix_cell = {gene_phenotype_score_matrix};
end

