
% parameter dist_metric 1:shortest path, 2:direct neighbor

function [gene_profile] = get_gene_profile(gene_id, g_p_network, SP, dist_metric)
len_phenotype = size(g_p_network, 2);
gene_profile = zeros(len_phenotype, 1);

for i = 1 : len_phenotype
    % finding disease phenotype genes.    
    sp_vec = SP(gene_id, g_p_network(:,i)>0);
    sp2 = -1*(sp_vec.*sp_vec);
    
    if dist_metric == 1
        topo_dis = sum(exp(sp2));
    elseif dist_metric == 2
        topo_dis = sum(exp(sp2(sp2 == -1)));         
    end    
    
    gene_profile(i, 1) = topo_dis;
end
