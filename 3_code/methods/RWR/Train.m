function [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train, matrix_cell_totrain)
    %cv_train_parameter_cell = {alpha; max_ites;indistct_num;length(indice_set)};
    %matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ncbi_gene_id;phenotype_id};
    %matrix_cell_totrain = MergeData(matrix_split_output_cell,i);
    parameter = cv_train_parameter_cell{1,1};
    lambda = parameter(1,1);
    max_ites = cv_train_parameter_cell{2,1};
    
    gene_phenotype_matrix_old = matrix_cell_totrain{1,1};
    phenotype_similarity_matrix = matrix_cell_train{2,1};
    ppi_matrix = matrix_cell_train{3,1};
    ncbi_gene_id = matrix_cell_train{4,1};
    phenotype_id = matrix_cell_train{5,1};
    
    W = Disisolate_ppi(ppi_matrix);
    D = zeros(size(W));
    for j = 1:size(W,1)
       D(j,j) = sum(W(j,:));
    end
    S = D^(-0.5) * W * D^(-0.5);


    array = 1./sum(gene_phenotype_matrix_old);
    array(isinf(array)) = 0;
    Y = gene_phenotype_matrix_old* diag(array);
    F_before = Y;
    for step = 1:20
      F_now = (1-lambda) * S * F_before + lambda * Y; 
      F_before = F_now;
    end
    gene_phenotype_score_matrix = F_now;


    learned_matrix_cell = {gene_phenotype_score_matrix};
end

