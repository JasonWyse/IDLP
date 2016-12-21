function [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train, matrix_cell_totrain)
    %cv_train_parameter_cell = {alpha; max_ites;indistct_num;length(indice_set)};
    %matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ncbi_gene_id;phenotype_id};
    %matrix_cell_totrain = MergeData(matrix_split_output_cell,i);
    beta = cv_train_parameter_cell{1,1}(1,1);

    max_ites = cv_train_parameter_cell{2,1};
    
    gene_phenotype_matrix_old = matrix_cell_totrain{1,1};
    phenotype_similarity_matrix = matrix_cell_train{2,1};
    ppi_matrix = matrix_cell_train{3,1};
    ncbi_gene_id = matrix_cell_train{4,1};
    phenotype_id = matrix_cell_train{5,1};
    
    gene_phenotype_score_matrix = zeros(size(gene_phenotype_matrix_old));
    array = sum(ppi_matrix,2);
    ppi_matrix_diag = diag(array);
    
    gene_phenotype_score_matrix = exp(-beta*(ppi_matrix_diag-ppi_matrix)) * gene_phenotype_matrix_old;
    
    learned_matrix_cell = {gene_phenotype_score_matrix};
end

