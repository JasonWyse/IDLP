% temp = Y_hat_folds_cell{1,6};
% [B,I] = sort(temp(:,1),'descend');
% [B1,I1] = sort(gene_phenotype_matrix_old(:,1),'descend');
% ans1 = sort(temp,'descend');
% ans2 = sort(gene_phenotype_matrix_old,'descend');
% [row,col,val] = find(ans1==1);
% [row1,col1,val1] = find(ans2 == 1);
% matrix_split_output_cell{1,1};

% fold_num = 5;
% statistics = zeros(5,5);
% for i=1:fold_num    
%     oneFold = matrix_split_output_cell{1,i};
%     otherFolds_cell = MergeData(matrix_split_output_cell,i);
%     otherFolds = otherFolds_cell{1,1};
%     oneFold_Y = Y_hat_folds_cell{1,i};
%     Y = 1./(1+exp((-15*oneFold_Y)+log(9999)));
%     oneFold_and_Y = oneFold.*(Y>0.9);
%     statistics(1,i) = nnz(oneFold);
%     statistics(2,i) = nnz(otherFolds);
%     statistics(3,i) = nnz(Y==1);
%     %tmp = oneFold_Y == 1;
%     %nnz(oneFold.*((Y>0.8)-otherFolds))
%     %isequal(otherFolds,otherFolds.*(oneFold_Y==1))
%     statistics(4,i) = nnz(otherFolds.*(Y==1));
%     statistics(5,i) = nnz(oneFold.*(Y>0.9));
%     statistics(6,i) = nnz(oneFold.*(Y==1));
% end
%[ Y_hat ] = Get_Y_hat_allFolds( gene_phenotype_matrix_old, phenotype_similarity_matrix );
Y_hat_allFolds = Y_hat_folds_cell{1,6};
Y = 1./(1+exp((-15*Y_hat_allFolds)+log(9999)));
%gene_phenotype_matrix_newAdded;
%gene_phenotype_matrix_old;
idx0 = sum(Y>0.98)>0;
nnz((Y>0.8).*gene_phenotype_matrix_newAdded)
nnz(gene_phenotype_matrix_old.*gene_phenotype_matrix_newAdded)
idx1 = sum(gene_phenotype_matrix_newAdded)>0;
idx2 = sum(gene_phenotype_matrix_old)>0;
nnz(idx1.*idx2)
singleton_disease = gene_phenotype_matrix_newAdded(:,sum(gene_phenotype_matrix_old)==0);
singleton_disease_Y = Y(:,sum(gene_phenotype_matrix_old)==0);
[B,I] = sort(singleton_disease_Y,'descend');
[r,c] = find(singleton_disease);
rank = zeros(1,length(c));
for i=1:length(c)
    [rank(i),~] = find(I(:,i)==r(i));
end










