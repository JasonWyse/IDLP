function [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train, matrix_cell_totrain)
    %cv_train_parameter_cell = {alpha; max_ites;indistct_num;length(indice_set)};
    %matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ncbi_gene_id;phenotype_id};
    %matrix_cell_totrain = MergeData(matrix_split_output_cell,i);
    alpha = cv_train_parameter_cell{1,1}(1,1);
    beta = cv_train_parameter_cell{1,1}(2,1);

    max_ites = cv_train_parameter_cell{2,1};
    
    gene_phenotype_matrix_old = matrix_cell_totrain{1,1};
    phenotype_similarity_matrix_norm = matrix_cell_train{7,1};
    ppi_matrix_norm = matrix_cell_train{8,1};
    
    
   [~,all_phenotype_num] = size(gene_phenotype_matrix_old);
    array = zeros(1,all_phenotype_num);
    array(1,1:end) = 1;
    p0 = diag(array);
 
    g0 = zeros(size(gene_phenotype_matrix_old));


    P_before = p0;
    G_before = g0;
    G_now = G_before;
    for i = 1:max_ites
        P_now = beta * phenotype_similarity_matrix_norm * P_before + (1-2*beta)* p0 + beta * gene_phenotype_matrix_old'* G_now; 
        G_now = alpha * ppi_matrix_norm * G_before + (1-2*alpha)*g0 + alpha*gene_phenotype_matrix_old*P_now;
        P_before = P_now;
        G_before = G_now;
    end
    
    gene_phenotype_score_matrix = G_now;
    learned_matrix_cell = {gene_phenotype_score_matrix};
end

