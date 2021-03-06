function [learned_matrix_cell] = Train2(cv_train_parameter_cell, matrix_cell_train, matrix_cell_totrain)
    %cv_train_parameter_cell = {alpha; max_ites;indistct_num;length(indice_set)};
    %matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ncbi_gene_id;phenotype_id};
    %matrix_cell_totrain = MergeData(matrix_split_output_cell,i);
    lambda = cv_train_parameter_cell{1,1}(1,1);

    max_ites = cv_train_parameter_cell{2,1};
    
    gene_phenotype_matrix_old = matrix_cell_totrain{1,1};
    phenotype_similarity_matrix_norm = matrix_cell_train{7,1};
    ppi_matrix_norm = matrix_cell_train{8,1};
    
    phenotype_gene_score_matrix = zeros(size(gene_phenotype_matrix_old,2),size(gene_phenotype_matrix_old,1));
    
    array = 1./sum(gene_phenotype_matrix_old);
    array(isinf(array)) = 0;
    gene_phenotype_matrix_old_norm = gene_phenotype_matrix_old* diag(array);
   
    R_before = phenotype_gene_score_matrix;
    for step = 1:max_ites
       R_now = lambda * phenotype_similarity_matrix_norm*R_before*ppi_matrix_norm...,
           + (1-lambda)*gene_phenotype_matrix_old_norm';
       R_before = R_now; 
    end
    gene_phenotype_score_matrix =  R_now';
    
    learned_matrix_cell = {gene_phenotype_score_matrix};
end

