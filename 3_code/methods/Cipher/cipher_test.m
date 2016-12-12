function cipher_test(dist_metric)

if strcmp(dist_metric,'sp') == 1
    dm = 1;
    fn = 'cipher_result_test_sp.mat';
elseif strcmp(dist_metric,'dn') == 1
    dm = 2;
    fn = 'cipher_result_test_dn.mat';
end

load ('DataSet.mat');
g_p_network = A_old;
phenotype_profile = phenotype_logistic;
[rows, cols] = size(g_p_network);

R = zeros(rows, cols);
R = CipherRank(phenotype_profile, g_p_network, 1, cols, ppi_sp, dm);

cROC = zeros(1,6);
[~,avgROC] = ROC_main(A_new',R');
cROC(1,:) = avgROC;

save (fn,'R','cROC');
end