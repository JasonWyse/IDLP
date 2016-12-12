rec = zeros(30,2);
load('Dataset.mat');
load('cipher_result_sp.mat');
rec(:,1) = recall(A_old',R');

load('cipher_result_test_sp.mat');
rec(:,2) = recall(A_new',R');

save('recall_sp.mat', 'rec');

load('cipher_result_dn.mat');
rec(:,1) = recall(A_old',R');

load('cipher_result_test_dn.mat');
rec(:,2) = recall(A_new',R');

save('recall_dn.mat', 'rec');