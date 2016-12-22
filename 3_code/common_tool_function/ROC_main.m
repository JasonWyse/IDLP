function [ROCn,avgROC] = ROC_main( phenotype_gene_matrix_groundTruth,phenotype_gene_matrix_prediction)
%R �� row_phenotype, column_gene
phenotype_gene_matrix_prediction(isnan(phenotype_gene_matrix_prediction)) = -1;
[rows, cols] = size(phenotype_gene_matrix_groundTruth);

ROCn = zeros(rows, 7);

[B, IX] = sort(phenotype_gene_matrix_prediction, 2, 'descend');%B�洢������������IX�洢����������,2:��������
clear R;
clear B;

topn = zeros(rows, cols);
for j = 1 : rows     %j: iterate phenotype 
    for k = 1 : cols%���ǰѰ�˳���ŵ�һ����Ԥ��Ե���һ����ǰ��һԽ��Խ��
        real_col = IX(j, k);
        if phenotype_gene_matrix_groundTruth(j, real_col) > 0
            topn(j, k) = 1;%
        else
            topn(j, k) = 0;
        end
    end

    ROCn(j,1) = ROC(topn(j,:), 20);
    ROCn(j,2) = ROC(topn(j,:), 50);
    ROCn(j,3) = ROC(topn(j,:), 100);
    ROCn(j,4) = ROC(topn(j,:), 300);
    ROCn(j,5) = ROC(topn(j,:), 500);
    ROCn(j,6) = ROC(topn(j,:), 1000);
    ROCn(j,7) = ROC(topn(j,:), cols);
end
%filtering out phenotypes that donot interact with any genes in new added
%gene-phenotype matrix;
avgROC = mean(ROCn(sum(ROCn,2)>0,:));

end
