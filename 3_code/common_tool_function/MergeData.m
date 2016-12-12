function [matrix_cell_train ] = MergeData( matrix_cell,index )
    matrix_cell_train = cell(size(matrix_cell,1),1);
    for i = 1: size(matrix_cell,1)
        matrix_temp = zeros(size(matrix_cell{i,1}));
        for j = 1:size(matrix_cell,2)
            if j~=index
               matrix_temp = matrix_temp + matrix_cell{i,j};
            end
        end
        matrix_cell_train{i,1} = matrix_temp;
    end
end

