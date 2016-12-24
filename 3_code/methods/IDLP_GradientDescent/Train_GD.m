function [learned_matrix_cell] = Train_GD(cv_train_parameter_cell, matrix_cell_train, ...,
                                    matrix_folds_cell_train,initialMatrix_cell,fold_idx)
   % cv_train_parameter_cell = {best_parameter_array; matrix_cv_split_idx; fold_num; 
   %                              distinct_parameter_num; evaluation_index_num};
   %matrix_cell_train = {gene_phenotype_matrix_old; phenotype_similarity_matrix; ppi_matrix; 
   %                     ncbi_gene_id; phenotype_id};
   %matrix_folds_cell_train = matrix_cell_train(matrix_cv_split_idx,1);
   %initialMatrix_cell = initialMatrix_cell = initialMatrice_cells{:,i};
   mu = cv_train_parameter_cell{1,1}(1,1);  
   nu = cv_train_parameter_cell{1,1}(2,1);
   zeta = cv_train_parameter_cell{1,1}(3,1);
   eta = cv_train_parameter_cell{1,1}(4,1);      
   max_ite = cv_train_parameter_cell{end-4+1,1};
   
   gene_phenotype_matrix_old = matrix_folds_cell_train{1,1};
   phenotype_similarity_matrix = matrix_cell_train{2,1};
   ppi_matrix = matrix_cell_train{3,1};
   Y_prince_folds_cell = matrix_cell_train{7,1};  
   
   W1 = Disisolate_ppi(ppi_matrix);
   D1 = diag(1./sqrt(sum(W1,2)));%    D1 = diag(sum(W1,2));
   S1_hat = D1*W1*D1;   %S1_hat = D1^(-0.5) * W1 * D1^(-0.5);
   W2 = phenotype_similarity_matrix; 
   D2 = diag(1./sqrt(sum(W2,2)));
   S2_hat = D2*W2*D2;  
   S1 = S1_hat;
   S2 = S2_hat;
   Y_init = initialMatrix_cell{1,1};
   Y = Y_init;   
   Y_hat = Y_prince_folds_cell{1,fold_idx}; 
  % Y_hat = 1./(1+exp((-15*Y_hat)+log(9999)));
   total_loss_vector = zeros(max_ite+1,1);
   total_loss_vector(1,1) =  Get_loss_function(Y,S1,S1_hat,S2,S2_hat,Y_hat,mu,nu,zeta,eta);
   for i=1:max_ite
       Y = Gradient_update_Y(Y,S1,S2,Y_hat,mu,zeta);
       S1 = Gradient_update_S1(Y,S1,S1_hat,nu);
       S2 = Gradient_update_S2(Y,S2,S2_hat,eta);
       total_loss_vector(i+1,1) =  Get_loss_function(Y,S1,S1_hat,S2,S2_hat,Y_hat,mu,nu,zeta,eta);
       
   end
   learned_matrix_cell = {Y;S1;S2};
end
function total_loss =  Get_loss_function(Y,S1,S1_hat,S2,S2_hat,Y_hat,mu,nu,zeta,eta)
    I1 = eye(size(S1));
    I2 = eye(size(S2));
    trace_loss = trace(Y'*(I1-S1)*Y) + trace(Y*(I2-S2)*Y');
    regularization_loss = (mu + zeta)*sum(sum(Y-Y_hat).^2) + nu*sum(sum(S1-S1_hat).^2)...,
        + eta*sum(sum(S2-S2_hat).^2);
    total_loss = trace_loss + regularization_loss;
    
end
function Y_new =  Gradient_update_Y(Y,S1,S2,Y_hat,mu,zeta)
    numerator = Y*S2 + S1*Y + (mu+zeta)*Y_hat;
    C = 2 + mu + zeta;
    denominator = C*Y;
    Y_new = Y.*sqrt(numerator./denominator);
    
end
function S1_new =  Gradient_update_S1(Y,S1,S1_hat,nu)
    numerator = 1/(2*nu)*(Y*Y') + S1_hat;    
    denominator = S1;
    S1_new = S1.*sqrt(numerator./denominator);
    
end
function S2_new =  Gradient_update_S2(Y,S2,S2_hat,eta)
    numerator = 1/(2*eta)*(Y'*Y) + S2_hat;    
    denominator = S2;
    S2_new = S2.*sqrt(numerator./denominator);    
end