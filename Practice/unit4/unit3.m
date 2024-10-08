%% submatrices
%A submatrix B of an mxn matrix is made of the elements of A such that two
%elements that belong to the same row or column of B also belong to the
%same row or column of A
clc
a = fix(rand(4,3)*10);
a(:,2);
a(3, :);
a(2:3, [1,3]);
b = a([2,4,1,3],[3,2,1]);

%% manipulating rows and columns
clc
%the matrix B=A(-[3,4],-[2,5]) is obtained from A by removing row 3 and 5
%and columns 2 and 5
A = fix(rand(6,7)*10);
[m,n] = size(A);
B = A(setdiff(1:m, [3,4]), setdiff(1:n, [2,5]));
B = A;
%B([3,4], [2,5]) = [];
C = A([3,4], [2,5])


%% Reshaping  a matrix
clc
a = reshape(1:30, 5,6);
b = a(setdiff(1:5, [3,4]), setdiff(1:6, [2,5]));
bb = a([1,2,5],[1,3,4,6]);
fliplr(a(1,:));
fliplr(a);
%common elements between flip(a(1,:) and 1:11 are 1, 6 and 7 11 and hence
%they are not displayed. The elements returned are 16 21 and 26 which are
%elements in flip(a(1,:) and not in 1:11
setdiff(flip(a(1,:)),1:11)
%setdiff always returns the results in ascending order hence [16 21 26]

%% Transpose and symmetry for real and compleex matrices
clc
%The transpose of an mxn matrix A is an mxn matrix b such tha aij=abji
%A symmetric matrix is a real matrix that is equal to its transpose

%SKEW-SYMMETRIC MATRIX,(antisymmetric), A=-A' meaning that after a mayrix
%is transposed, it becomes a negative of itself
%all the diagonal elements must be zero, a square matrix and the elemnts
%across the diagonal are negative of each other.

A = [0 2 -4; -2 0 6; 4 -6 0];
A';
isequal(A', -A);
%SKEW-HERMITIAN MATRIX, a square complex matrix A that Aconjugate=-A. if
%you take the transpose of A and conjugate each element(replace every
%element with its complex conjugete), the result is the negative of the
%original matrix
Ah = [1+j 1; -2j 2+3j];
isequal(Ah', -conj(Ah));
j=sqrt(-1);
a = [1+j, 1; -2*j, (2+3*j)];

A = [ 2-1i 3+2i; -2-1i 0 4i; -3-2i -4i -1i];  % Example skew-Hermitian matrix

if isequal(A', -conj(A))
    disp('The matrix is skew-Hermitian.');
else
    disp('The matrix is not skew-Hermitian.');
end

%% STORAGE AND SPARCITY
%by default, matlab handles matrices as dense, each elemnt is stored which
%can lead to ghigh memory usage for large matrices hence why sparse
%matrices are important. for matrices that contain a significnat number of
%zeros, one can use sparse storage whee only non-zero elements are stored
clc

rows = [1 3 4];
cols = [1 3 2];
values = [10 20 30];
A = sparse(rows, cols, values, 5, 5);
B = sparse(A);
%Checking for sparcity of a matrix using nnz(number of non-zero elements
density = nnz(A) / numel(A);   % Density of non-zero elements

%Create a large matrix
%define the dimensions of the sparce matrix
n = 1000;   
% randi(n, 10000, 1): This generates a column vector of 10,000 random integers between 1 and 1000 (inclusive). This vector represents the row indices of the non-zero elements in the sparse matrix.
% randi(n, 10000, 1): Similarly, this generates another column vector of 10,000 random integers between 1 and 1000, representing the column indices of the non-zero elements.
% rand(10000, 1): This generates a column vector of 10,000 random numbers uniformly distributed between 0 and 1. These numbers are the values of the non-zero elements in the sparse matrix.

% S = sparce(randi(n, 10000, 1), randi(n, 10000, 1), rand(10000, ) n, n)
S = sparse(randi(n, 10000, 1), randi(n, 10000, 1), rand(10000, 1), n, n)