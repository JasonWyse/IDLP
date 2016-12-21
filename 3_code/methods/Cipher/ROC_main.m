function [ROCn,avgROC] = ROC_main( HPG,R)
%R ： row_phenotype, column_gene
R(isnan(R)) = -1;
[rows, cols] = size(HPG);

ROCn = zeros(rows, 6);

[B, IX] = sort(R, 2, 'descend');%B存储的是排序结果，IX存储的是索引列
clear R;
clear B;

topn = zeros(rows, cols);
for j = 1 : rows     %j: iterate phenotype 
    for k = 1 : cols%就是把按顺序排的一行中预测对的置一，靠前的一越多越好
	real_col = IX(j, k);
        if (HPG(j, real_col) > 0)
            topn(j, k) = 1;%
        else
            topn(j, k) = 0;
        end
    end

    ROCn(j,1) = ROC(topn(j,:), 50);
    ROCn(j,2) = ROC(topn(j,:), 100);
    ROCn(j,3) = ROC(topn(j,:), 300);
    ROCn(j,4) = ROC(topn(j,:), 500);
    ROCn(j,5) = ROC(topn(j,:), 1000);
    ROCn(j,6) = ROC(topn(j,:), cols);
end
avgROC = mean(ROCn);

% sample_size = 50;
% distri = zeros(sample_size, 1+6);
% for i = 1 : sample_size
%     distri(i, 1) = i / sample_size;
% end
% 
% for i = 1 : sample_size
%     distri(i, 2) = sum(ROCn(:,1)>=distri(i, 1));
%     distri(i, 3) = sum(ROCn(:,2)>=distri(i, 1));
%     distri(i, 4) = sum(ROCn(:,3)>=distri(i, 1));
%     distri(i, 5) = sum(ROCn(:,4)>=distri(i, 1));
%     distri(i, 6) = sum(ROCn(:,5)>=distri(i, 1));
%     distri(i, 7) = sum(ROCn(:,6)>=distri(i, 1));
% end
end