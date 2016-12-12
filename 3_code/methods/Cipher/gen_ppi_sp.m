load 'DataSet.mat';
ppi_sp = graphallshortestpaths(sparse(ppi_network));
save('DataSet.mat','ppi_sp','A_old','A_new','phenotype_logistic','ppi_network');