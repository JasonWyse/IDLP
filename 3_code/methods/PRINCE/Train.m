function [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train, matrix_cell_totrain)
    %cv_train_parameter_cell = {alpha; max_ites;indistct_num;length(indice_set)};
    %matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ncbi_gene_id;phenotype_id};
    %matrix_cell_totrain = MergeData(matrix_split_output_cell,i);
    parameter = cv_train_parameter_cell{1,1};
    alpha = parameter(1,1);
    max_ites = cv_train_parameter_cell{2,1};
    
    gene_phenotype_matrix_old = matrix_cell_totrain{1,1};
    phenotype_similarity_matrix = matrix_cell_train{2,1};
    ppi_matrix = matrix_cell_train{3,1};
    ncbi_gene_id = matrix_cell_train{4,1};
    phenotype_id = matrix_cell_train{5,1};
    
    W = Disisolate_ppi(ppi_matrix);
    D = diag(sum(W,2));
    S = D^(-0.5) * W * D^(-0.5);

    all_gene_num = length(ncbi_gene_id);
    all_phenotype_num = length(phenotype_id);
    Y_temp = zeros(length(ncbi_gene_id),length(phenotype_id));
    for i = 1 : all_gene_num
        tmp = repmat(gene_phenotype_matrix_old(i,:),all_phenotype_num,1);
        temp2 = tmp.*phenotype_similarity_matrix;
        Y_temp(i,:) = max(temp2,[],2)';

    end
    Y = 1./(1+exp((-15*Y_temp)+log(9999)));
    
    F_before = Y;
    for step = 1:max_ites
      F_now = alpha * S * F_before + (1-alpha) * Y; 
      F_before = F_now;
    end
    gene_phenotype_score_matrix = F_now;
    learned_matrix_cell = {gene_phenotype_score_matrix};
end

