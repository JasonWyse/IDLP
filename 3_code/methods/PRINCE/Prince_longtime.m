clear;
%refer to : Associating genes and Protein Complexes with Disease via Network Progagation
path(path,'../../../2_useful_data');
path(path,'../../common_tool_function');
GP_file_name = 'G_P_network_mappingkey13_2015_8&2016_12.mat';

alpha = 0.5;
load(GP_file_name,'gene_phenotype_matrix_old','gene_phenotype_matrix_newAdded',...,
    'phenotype_similarity_matrix','ncbi_gene_id','phenotype_id','ppi_matrix');
   
W = Disisolate_ppi(ppi_matrix);
D = diag(sum(W,2));
S = D^(-0.5) * W * D^(-0.5);
gene_phenotype_score_matrix = zeros(size(gene_phenotype_matrix_old));
for i = 1:length(phenotype_id)
    Y = zeros(length(ncbi_gene_id),1);
    for j = 1:length(ncbi_gene_id)
        [phenotypes] = find(gene_phenotype_matrix_old(j,:)==1);
        if(isempty(phenotypes))
             Y(j,1) = 0;
        else
            Y(j,1) = max(phenotype_similarity_matrix(i,phenotypes));
        end
    end
    Y = 1./(1+exp((-15*Y)+log(9999)));
    F_before = Y;
    for step = 1:20
       F_now = alpha * S * F_before + (1-alpha) * Y;
       F_before = F_now;
    end
    gene_phenotype_score_matrix(:,1) = F_now;
end

gene_phenotype_score_matrix%save('gene_phenotype_score_matrix.mat','gene_phenotype_score_matrix');
[~,cRoc] = ROC_main(gene_phenotype_matrix_newAdded',gene_phenotype_score_matrix');

rmpath('../../../2_useful_data');
rmpath('../../common_tool_function');
