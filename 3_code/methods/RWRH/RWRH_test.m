function [] = RWRH_test(phenotype_gene_new,phenotype_gene_score_matrix )
    path(path,'../../common_tool_function');
    cROC = zeros(1,6);
    [~,avgROC] = ROC_main(phenotype_gene_new,phenotype_gene_score_matrix);
    cROC(1,:) = avgROC;

    save ('result.mat','phenotype_gene_score_matrix','cROC');
end

