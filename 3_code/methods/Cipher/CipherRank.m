function [R] = CipherRank(phenotype_profile, g_p_network, read_begin, read_end, SP, dist_metric)
    R = zeros(size(g_p_network,1), read_end - read_begin + 1);

    [allGene_num, ~] = size(g_p_network);
    allgene_profile = zeros(size(g_p_network'));
    for i = 1:allGene_num 
        allgene_profile(:,i) = get_gene_profile(i, g_p_network, SP, dist_metric);
    end 
    matrix_conbine = [allgene_profile,phenotype_profile(:,read_begin:read_end)];
    pearson_matrix = corrcoef(matrix_conbine);
    R = pearson_matrix(1:allGene_num,allGene_num+1:end); 
end
