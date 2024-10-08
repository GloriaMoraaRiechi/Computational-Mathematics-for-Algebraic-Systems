
%distributes or updates values to specific locations in a sparse matrix or
%data structure

%initialize a sparse matrix
A = sparse(4, 4);

%indices where updates should be scattered
rowIndices = [1, 2, 3, 4];
colIndices = [1, 2, 4, 5];
values = [10, 20, 30, 40];

%ScatterU operation
for i = 1:length(values)
    A(rowIndices(i), colIndices(i)) = values(i);
end

full(A)
