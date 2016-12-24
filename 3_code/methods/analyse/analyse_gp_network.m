clear;
%refer to : Genome-wide inferring gene¨Cphenotype relationship by walking on the heterogeneous network.pdf
path(path,'../../../2_useful_data');
path(path,'../../common_tool_function');
file_date_time = '2015_8&2016_12';
%ppi_key_word = 'changeppi_';
ppi_key_word = '';

GP_file_name = ['G_P_network_' ppi_key_word 'mappingkey13_' file_date_time '.mat'];

load(GP_file_name,'gene_phenotype_matrix_old', 'gene_phenotype_matrix_newAdded', 'phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','phenotype_id','ppi_matrix','gene_phenotype_matrix_new');
[~,phenotypes] = find(gene_phenotype_matrix_newAdded == 1);
phenotypes_unique = unique(phenotypes);

array_phenotype_totally_new = zeros(1,1);
array_phenotype_exist_old = zeros(1,1);
size1 = 0;
size2 = 0;
for i = 1:length(phenotypes_unique)
   if sum(gene_phenotype_matrix_old(:,phenotypes_unique(i))) == 0
       array_phenotype_totally_new(size1+1) = phenotypes_unique(i);
       size1 = size1 + 1;
   else
       array_phenotype_exist_old(size2+1) = phenotypes_unique(i);
       size2 = size2 + 1;
   end
end

array_phenotype_exist_old_genes = sum(gene_phenotype_matrix_old(:,array_phenotype_exist_old));
%[array_phenotype_exist_old_genes_sort,~] = sort(array_phenotype_exist_old_genes);

array_phenotype_exist_old_genes(2,:) = sum(gene_phenotype_matrix_newAdded(:,array_phenotype_exist_old));
%[array_phenotype_exist_old_genes_sort(2,:),~] = sort(array_phenotype_exist_old_genes);

%exist_old
gene_phenotype_matrix_newAdded_temp1 = zeros(size(gene_phenotype_matrix_newAdded));
gene_phenotype_matrix_newAdded_temp1(:,array_phenotype_exist_old) = gene_phenotype_matrix_newAdded(:,array_phenotype_exist_old);
%total_new
gene_phenotype_matrix_newAdded_temp2 = zeros(size(gene_phenotype_matrix_newAdded));
gene_phenotype_matrix_newAdded_temp2(:,array_phenotype_totally_new) = gene_phenotype_matrix_newAdded(:,array_phenotype_totally_new);

ppi_key_word = 'exist_old_';
GP_file_name = ['G_P_network_' ppi_key_word 'mappingkey13_' file_date_time '.mat'];
gene_phenotype_matrix_newAdded = gene_phenotype_matrix_newAdded_temp1;
save(GP_file_name,'gene_phenotype_matrix_old', 'gene_phenotype_matrix_newAdded', 'phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','phenotype_id','ppi_matrix','gene_phenotype_matrix_new');

ppi_key_word = 'total_new_';
GP_file_name = ['G_P_network_' ppi_key_word 'mappingkey13_' file_date_time '.mat'];
gene_phenotype_matrix_newAdded = gene_phenotype_matrix_newAdded_temp2;
save(GP_file_name,'gene_phenotype_matrix_old', 'gene_phenotype_matrix_newAdded', 'phenotype_similarity_matrix'...,
    ,'ncbi_gene_id','phenotype_id','ppi_matrix','gene_phenotype_matrix_new');