clear;
%refer to : Walking the interactome for Prioritization of Candidate Disease Genes
path(path,'../../../2_useful_data');
path(path,'../../common_tool_function');
GP_file_name = 'G_P_network_mappingkey13_2015_8&2016_12.mat';

lambda1 = 0.5;

load(GP_file_name,'gene_phenotype_matrix_old','gene_phenotype_matrix_newAdded','phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','phenotype_id','ppi_matrix');
W = Disisolate_ppi(ppi_matrix);
D = zeros(size(W));
for j = 1:size(W,1)
   D(j,j) = sum(W(j,:));
end
S = D^(-0.5) * W * D^(-0.5);



Y = gene_phenotype_matrix_old* diag(1./sum(gene_phenotype_matrix_old));
F_before = Y;
for step = 1:20
  F_now = (1-lambda1) * S * F_before + lambda1 * Y; 
  F_before = F_now;
end
gene_phenotype_score_matrix = F_now;

save('gene_phenotype_score_matrix.mat','gene_phenotype_score_matrix');
RWR_test(gene_phenotype_matrix_newAdded',gene_phenotype_score_matrix');

rmpath('../../../2_useful_data');
rmpath('../../common_tool_function');
