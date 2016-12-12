function cipher_train(dist_metric)

if strcmp(dist_metric,'sp') == 1
    dm = 1;
    fn = 'cipher_result_sp.mat';
elseif strcmp(dist_metric,'dn') == 1
    dm = 2;
    fn = 'cipher_result_dn.mat';
end

load ('DataSet.mat');
g_p_network = A_old;
phenotype_profile = phenotype_logistic;
[rows, cols] = size(g_p_network);

R = zeros(rows, cols);
loop=5;
for i = 1 : loop
    disp(['Fold:' num2str(i) 'th is running... ']);
    tic;
    read_begin = round((i - 1) * cols / loop) + 1;
    read_end = round(i * cols / loop);
    
    tmp_buffer = g_p_network(:, read_begin:read_end);
    g_p_network(:,read_begin:read_end) = 0;
    
    Rslt = CipherRank(phenotype_profile, g_p_network, read_begin, read_end, ppi_sp, dm);
    
    R(:,read_begin:read_end) = Rslt;
    g_p_network(:,read_begin:read_end) = tmp_buffer;
    toc;
end

cROC = zeros(2,6);
[~,avgROC] = ROC_main(A_old',R');
cROC(1,:) = avgROC;
[~,avgROC] = ROC_main(A_new',R');
cROC(2,:) = avgROC;

save (fn,'R','cROC');
end