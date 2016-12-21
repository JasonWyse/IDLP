function [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train, matrix_cell_totrain)
    %cv_train_parameter_cell = {method_str; max_ites;length(indice_set)};
    %matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ppi_sp;ncbi_gene_id;phenotype_id};
    %matrix_cell_totrain = MergeData(matrix_split_output_cell,i);
    method_str = cv_train_parameter_cell{1,1};
    max_ites = cv_train_parameter_cell{2,1};
    
    gene_phenotype_matrix_old = matrix_cell_totrain{1,1};
    phenotype_similarity_matrix = matrix_cell_train{2,1};
    ppi_matrix = matrix_cell_train{3,1};
    ppi_sp = matrix_cell_train{4,1};
    ncbi_gene_id = matrix_cell_train{5,1};
    phenotype_id = matrix_cell_train{6,1};
    
    if strcmp(method_str,'sp') == 1
        dm = 1;
        fn = 'cipher_result_sp.mat';
    elseif strcmp(method_str,'dn') == 1
        dm = 2;
        fn = 'cipher_result_dn.mat';
    end

    [rows, cols] = size(gene_phenotype_matrix_old);
    gene_phenotype_score_matrix = zeros(rows, cols);
    for i = 1 : max_ites
        disp(['Fold:' num2str(i) 'th is running... ']);
        tic;
        read_begin = round((i - 1) * cols / loop) + 1;
        read_end = round(i * cols / loop);

        tmp_buffer = gene_phenotype_matrix_old(:, read_begin:read_end);
        gene_phenotype_matrix_old(:,read_begin:read_end) = 0;

        Rslt = CipherRank(phenotype_similarity_matrix, gene_phenotype_matrix_old, read_begin, read_end, ppi_sp, dm);

        gene_phenotype_score_matrix(:,read_begin:read_end) = Rslt;
        g_p_network(:,read_begin:read_end) = tmp_buffer;
        toc;
    end

cROC = zeros(2,6);
[~,avgROC] = ROC_main(A_old',gene_phenotype_score_matrix');
cROC(1,:) = avgROC;
[~,avgROC] = ROC_main(A_new',gene_phenotype_score_matrix');
cROC(2,:) = avgROC;

save (fn,'R','cROC');

    learned_matrix_cell = {gene_phenotype_score_matrix};
end

