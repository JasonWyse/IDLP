clear;
%refer to : Associating genes and Protein Complexes with Disease via Network Progagation
path(path,'../../../2_useful_data');
path(path,'../../common_tool_function');
GP_file_name = 'G_P_network_mappingkey13_2015_8&2016_12.mat';


load(GP_file_name,'gene_phenotype_matrix_old','gene_phenotype_matrix_newAdded',...,
    'phenotype_similarity_matrix','ncbi_gene_id','phenotype_id','ppi_matrix');
   
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
save('Y.mat','Y_temp');
load('Y.mat','Y_temp');
%-----load is equal to above lines

%Y = 1./(1+exp((-15*Y_temp)+log(9999)));
alpha = 0.9;
Y = Y_temp;
F_before = Y;
for step = 1:20
  F_now = alpha * S * F_before + (1-alpha) * Y; 
  F_before = F_now;
end
gene_phenotype_score_matrix = F_now;%save('gene_phenotype_score_matrix.mat','gene_phenotype_score_matrix');
[~,cRoc] = ROC_main(gene_phenotype_matrix_newAdded',gene_phenotype_score_matrix');

save ('result.mat','phenotype_gene_score_matrix','cROC');
%Prince_test(gene_phenotype_matrix_newAdded',gene_phenotype_score_matrix');
rmpath('../../../2_useful_data');
rmpath('../../common_tool_function');
