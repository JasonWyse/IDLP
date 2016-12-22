function [ Y_hat ] = Initialize_Y( gene_phenotype_matrix)
   array = 1./sum(gene_phenotype_matrix);
   array(isinf(array)) = 0;
   Y = gene_phenotype_matrix* diag(array);
    
   Y_hat = Y; 

end

