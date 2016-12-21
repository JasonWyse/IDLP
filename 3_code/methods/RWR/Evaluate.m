function [ AUC_vec ] = Evaluate( learned_matrix_cell, matrix_validation_cell,index_num )

    gene_phenotype_score_matrix = learned_matrix_cell{1,1};
    %matrix_cell_test = {gene_phenotype_matrix_newAdded};
    gene_phenotype_matrix_newAdded =  matrix_validation_cell{1,1};
   % path(path,'../../common_tool_function');
    cROC = zeros(1,index_num);
    phenotype_gene_newAdded =  gene_phenotype_matrix_newAdded';
    phenotype_gene_score_matrix = gene_phenotype_score_matrix';
    [~,avgROC] = ROC_main(phenotype_gene_newAdded, phenotype_gene_score_matrix);
    cROC(1,:) = avgROC;

    %save ('result.mat','phenotype_gene_score_matrix','cROC');
    AUC_vec =cROC';

end

