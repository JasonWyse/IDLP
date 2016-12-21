function [learned_matrix_cell] = Train(cv_train_parameter_cell, matrix_cell_train, matrix_cell_totrain)
    %cv_train_parameter_cell = {alpha; max_ites;indistct_num;length(indice_set)};
    %matrix_cell_train = {gene_phenotype_matrix_old;phenotype_similarity_matrix;ppi_matrix;ncbi_gene_id;phenotype_id};
    %matrix_cell_totrain = MergeData(matrix_split_output_cell,i);
    lambda = cv_train_parameter_cell{1,1}(1,1);
    eta = cv_train_parameter_cell{1,1}(2,1);
    gamma = cv_train_parameter_cell{1,1}(3,1);

    max_ites = cv_train_parameter_cell{2,1};
    
    gene_phenotype_matrix_old = matrix_cell_totrain{1,1};
    phenotype_similarity_matrix = matrix_cell_train{2,1};
    ppi_matrix = matrix_cell_train{3,1};
    ncbi_gene_id = matrix_cell_train{4,1};
    phenotype_id = matrix_cell_train{5,1};
    
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
    
    array = 1./sum(gene_phenotype_matrix_old);
    array(isinf(array)) = 0;
    U0_matrix = gene_phenotype_matrix_old* diag(array);
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
    
    learned_matrix_cell = {gene_phenotype_score_matrix};
end

