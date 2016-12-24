clear;
path(path,'../../../2_useful_data');
path(path,'../../common_tool_function');
file_date_time = '2015_8&2016_12';

GP_file_name = ['G_P_network_mappingkey13_' file_date_time '.mat'];

load(GP_file_name);
%%%%%%%%%%%%%%%%%%%%parameter%%%%%%%%%%%%%%%%%%%%
rate1 = 0.1;    %remove known ppi
rate2 = 0.1;    %add ppi
%%%%%%%%%%%%%%%%%%%%parameter%%%%%%%%%%%%%%%%%%%%
all_ppi_num = nnz(ppi_matrix);
remove_num = floor(all_ppi_num * rate1*0.5);

[r,c] = find(ppi_matrix == 1);
array = randperm(all_ppi_num);
for i = 1:remove_num
   row = r(array(i)) ;
   col = c(array(i));
   ppi_matrix(row,col) = 0;
   ppi_matrix(col,row) = 0;

end

array1 = randperm(size(ppi_matrix,1));
array2 = randperm(size(ppi_matrix,1));
for i = 1:remove_num
   row = array1(i) ;
   col = array2(i);
   ppi_matrix(row,col) = 1;
   ppi_matrix(col,row) = 1;

end
GP_file_name = ['G_P_network_changeppi_mappingkey13_' file_date_time '.mat'];


save(GP_file_name,'gene_phenotype_matrix_new','gene_phenotype_matrix_newAdded'...,
    ,'gene_phenotype_matrix_old','ncbi_gene_id','phenotype_id',...,
    'phenotype_similarity_matrix','ppi_matrix');