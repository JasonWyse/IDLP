function [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train, matrix_cell_totrain)
    %cv_train_parameter_cell = {alpha; max_ites;indistct_num;length(indice_set)};
    %matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ncbi_gene_id;phenotype_id};
    %matrix_cell_totrain = MergeData(matrix_split_output_cell,i);
    beta = cv_train_parameter_cell{1,1}(1,1);
    
    gene_phenotype_matrix_old = matrix_cell_totrain{1,1};
    ppi_matrix_minus = matrix_cell_train{7,1};
    
    gene_phenotype_score_matrix = zeros(size(gene_phenotype_matrix_old));
    
    gene_phenotype_score_matrix = exp(-beta*(ppi_matrix_minus)) * gene_phenotype_matrix_old;
    
    learned_matrix_cell = {gene_phenotype_score_matrix};
end

