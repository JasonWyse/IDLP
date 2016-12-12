function [rec] = recall( HPG,R)

R(isnan(R)) = -1;
[rows, cols] = size(HPG);
rec = zeros(30, 1);

[~, IX] = sort(R, 2, 'descend');%B�洢������������IX�洢����������
clear R;

topn = zeros(rows, cols);
for j = 1 : rows
    for k = 1 : cols%���ǰѰ�˳���ŵ�һ����Ԥ��Ե���һ����ǰ��һԽ��Խ��
	real_col = IX(j, k);
        if (HPG(j, real_col) > 0)
            topn(j, k) = 1;%
        else
            topn(j, k) = 0;
        end
    end
end
total = nnz(topn);
for i = 10:10:300
    tmp = topn(:,1:i);
    rec(i/10,1) = nnz(tmp)/total;
end
end