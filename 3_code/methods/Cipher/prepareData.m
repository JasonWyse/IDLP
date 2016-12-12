clear;
path(path,'../../../2_useful_data');
load('G_PInfo_mapping13_2015_8.mat','gene_phenotype_matrix','phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','omim_gene_id','phenotype_id','ppi_ppi');

gene_phenotype_matrix1 = gene_phenotype_matrix;
ncbi_gene_id1 = ncbi_gene_id;
phenotype_id1 = phenotype_id;
load('G_PInfo_mapping13_2016_12_2.mat','gene_phenotype_matrix','phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','omim_gene_id','phenotype_id','ppi_ppi');
gene_phenotype_matrix2 = gene_phenotype_matrix;
ncbi_gene_id2 = ncbi_gene_id;
phenotype_id2 = phenotype_id;

phenotype_id = unique([phenotype_id1;phenotype_id2]);
ncbi_gene_id = unique([ncbi_gene_id1;ncbi_gene_id2]);

A_old = zeros(length(ncbi_gene_id),length(phenotype_id));
A_new = zeros(length(ncbi_gene_id),length(phenotype_id));

[r,c] = find(gene_phenotype_matrix1 == 1);
for i = 1:length(r)
   gene =  ncbi_gene_id1(r(i));
   phenotype = phenotype_id1(c(i));
   [gene_index,~] = find(ncbi_gene_id == gene);
   [phenotype_index,~] = find(phenotype_id == phenotype);
   A_old(gene_index,phenotype_index) = 1;
end

[r,c] = find(gene_phenotype_matrix2 == 1);
for i = 1:length(r)
   gene =  ncbi_gene_id2(r(i));
   phenotype = phenotype_id2(c(i));
   [gene_index,~] = find(ncbi_gene_id == gene);
   [phenotype_index,~] = find(phenotype_id == phenotype);
   A_new(gene_index,phenotype_index) = 1;
end

path(path,'../../../1_prepare_data/ppi');
load('ppi_ppi_2015_8.mat','ppi_index','ppi_ppi');
ppi_network = size(length(ncbi_gene_id),length(ncbi_gene_id));
[~,ia,ib] = intersect(ppi_index,ncbi_gene_id);
ppi_network(ib,ib) = ppi_ppi(ia,ia);

phenotype_id_origin = phenotype_id;
load('phenotype_similarity_2015_12.mat','phenotype_id','phenotype_similarity_matrix');
phenotype_logistic = size(length(phenotype_id_origin),length(phenotype_id_origin));
[~,ia,ib] = intersect(phenotype_id,phenotype_id_origin);
phenotype_logistic(ib,ib) = phenotype_similarity_matrix(ia,ia);

save('DataSet.mat','A_new','A_old','phenotype_logistic','ppi_network');

gene_phenotype_matrix_old = A_old;
gene_phenotype_matrix_new = A_new;
phenotype_id = phenotype_id_origin;
phenotype_similarity_matrix = phenotype_logistic;
ppi_matrix = ppi_network;
save('G_P_network_mappingkey13_2015_8&2016_12.mat','gene_phenotype_matrix_old','gene_phenotype_matrix_new',...,
    'ncbi_gene_id','phenotype_id','phenotype_similarity_matrix','ppi_matrix');


rmpath('../../../2_useful_data');
rmpath('../../../1_prepare_data/ppi');
