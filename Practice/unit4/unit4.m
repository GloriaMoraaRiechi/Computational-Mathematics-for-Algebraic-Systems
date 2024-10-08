%% VECTOR AND MATRIX ADDITION AND MULTIPLICATION
%% matrix scatter operation

clc
A = zeros(5, 5);

values = [10, 20, 30];    %values that need to be scattered
row_index = [ 1, 3, 5];   %Indices where values need to be scattered
col_index =  [2, 4, 1];

% Scatter operation: assigning values to the specified row and column indices
for i = 1:length(values)
    A(row_index(i), col_index(i)) = values(i);
end

disp('matrix after scatter operation: ');
disp(A)

%% matrix gather operation
clc
A = [10 20 30 40; 50 60 70 80; 90 100 110 120; 130 140 150 160];

%Rows and columns to gather values from
row_index = [ 1, 3, 4];
col_index = [2, 4, 1];

%Preallocate for the gathered values
gathered_values(i) = A(row_index(i), col_index(i));

for i=1:length(row_index)
    gathered_values(i) = A(row_index(i), col_index(i));
end

disp('gathered values: ');
disp(gathered_values);

%gatherU Function
%takes a full-storage vector u as input and returns a sparse
%representation of the nonzero elements of u. The output w is a two column
%matrix where the fisrt column stores indices of the non-zero elements and
%the second column stores the corresponding values of those elements.
u = [0, 5, 0, 0, 3, 0, 7];

function w = GatherU (u)
% u is full-storage and w is "sparse format" with two columns (index, value)
% the index of the nonzero element and the value
n = length(u) ; nnz = 0 ;
for i = 1 : n % loop over all elements in u
    if u(i) ~= 0 % ok if u is made of integers.
% If u is a floating point number, then must use abs(u(i))>10*eps
        nnz = nnz + 1 ;
        w(nnz,1) = i ;
        w(nnz,2) = u(i) ;
    end
end
end
%% 
clc
