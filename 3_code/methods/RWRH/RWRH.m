clear;
%refer to : Genome-wide inferring gene¨Cphenotype relationship by walking on the heterogeneous network.pdf
path(path,'../../../2_useful_data');
GP_file_name = 'G_P_network_mappingkey13_2015_8&2016_12.mat';

lambda = 0.5;
eta = 0.5;
gamma = 0.7;

load(GP_file_name,'gene_phenotype_matrix_old','gene_phenotype_matrix_newAdded',...,
    'phenotype_similarity_matrix','ncbi_gene_id','ppi_matrix');

adjacency_matrix = [ppi_matrix,gene_phenotype_matrix_old;gene_phenotype_matrix_old',phenotype_similarity_matrix];

M_GP = zeros(size(gene_phenotype_matrix_old));
temp_matrix = zeros(size(gene_phenotype_matrix_old));
for i = 1:size(gene_phenotype_matrix_old,1)
   temp_matrix(i,:) = sum(gene_phenotype_matrix_old(i,:)); 
end
M_GP = lambda * gene_phenotype_matrix_old./temp_matrix;
M_GP(isnan(M_GP)) = 0;


M_PG = zeros(size(gene_phenotype_matrix_old'));
temp_matrix = zeros(size(gene_phenotype_matrix_old'));
for i = 1:size(gene_phenotype_matrix_old,2)
   temp_matrix(i,:) = sum(gene_phenotype_matrix_old(:,i)); 
end
M_PG = lambda * gene_phenotype_matrix_old'./temp_matrix;
M_PG(isnan(M_PG)) = 0;

M_G = zeros(size(gene_phenotype_matrix_old,1),size(gene_phenotype_matrix_old,1));
for i = 1:size(gene_phenotype_matrix_old,1)
   sumB =  sum(gene_phenotype_matrix_old(i,:));
   sumA =  sum(ppi_matrix(i,:));
   if(sumB == 0)
       M_G(i,:) = ppi_matrix(i,:)./sumA;
   else
       M_G(i,:) = (1-lambda) * ppi_matrix(i,:)./sumA;
   end
    
end
M_G(isnan(M_G)) = 0;

M_P = zeros(size(gene_phenotype_matrix_old,2),size(gene_phenotype_matrix_old,2));
for i = 1:size(gene_phenotype_matrix_old,2)
   sumB =  sum(gene_phenotype_matrix_old(:,i));
   sumA =  sum(phenotype_similarity_matrix(i,:));
   if(sumB == 0)
       M_P(i,:) = phenotype_similarity_matrix(i,:)./sumA;
   else
       M_P(i,:) = (1-lambda) * phenotype_similarity_matrix(i,:)./sumA;
   end
    
end
M_P(isnan(M_P)) = 0;

M = [M_G,M_GP;M_PG,M_P];
[B,I] = sort(phenotype_similarity_matrix,1,'descend');

U0_matrix = gene_phenotype_matrix_old* diag(1./sum(gene_phenotype_matrix_old));
V0_matrix = zeros(size(M_P));
tmpI = (I(2:6,:));
[r,c] = size(tmpI);
col1 = reshape(tmpI,r*c,1);
tmp2 = repmat(I(1,:),5,1);
[r,c] = size(tmp2);
col2 = reshape(tmp2,r*c,1);
V0_matrix(sub2ind(size(V0_matrix),col1,col2)) = 0.2;
P0_matrix = [(1-eta)*U0_matrix;eta* V0_matrix];
P_before = P0_matrix;
    
for step = 1 : 30
   P_now = (1-gamma)*M'*P_before + gamma*P0_matrix;   
   P_before = P_now;
end
allGene_num = length(ncbi_gene_id);
gene_phenotype_score_matrix = P_now(1:allGene_num,:);

save('gene_phenotype_score_matrix.mat','gene_phenotype_score_matrix');

RWRH_test(gene_phenotype_matrix_newAdded',gene_phenotype_score_matrix');

rmpath('../../../2_useful_data');