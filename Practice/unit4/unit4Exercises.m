clc
A1 = reshape(1:16, 4,4)'
%getting a 2x2 top left submatrix of A1
B1 = A1(1:2, 1:2);

B2 = A1(:, [2,3,2,1,4]);

B3 = diag(diag(A1(1:3, 2:4)))

%% ABC=A(BC)
clc
A = randi(10,3,4);
B = randi(10,4,5);
C = randi(10,5,7);
prod = A * B * C ;     %results to a 3x7 matrix
prod2 = A * (B * C);

ans1 = (A*B) * 3;
ans2 = A * (B*3);

%%creating a sparse matrix
rows = 4;
columns = 4;
density = 0.3;
A3 = sprand(rows, columns, density);
size(A3);
full(A3);

B3 = sprand(rows, columns, density);
full(B3);


A3 * B3;

%% 
clc
%Define the matrix dimensions, number of rows and columns
n = 10;

%Define the densities
densityA = 0.1;
densityB = 0.05;

%Generate the random matrices
A = sprand(n, n, densityA);
full(A);

B = sprand(n, n, densityB);
full(B);

Cs = A + B
%count the number of non-zero elements in the sum
numOfnonzeroCs = nnz(Cs);
Cp = A * B
numOfnonzeroCp = nnz(Cp);
%Find the total number of elements in the matrix
totalElements = n * n;

%calculate the density of the sum and product
densityCs = numOfnonzeroCs / totalElements;
densityCp = numOfnonzeroCp / totalElements;
fprintf('The density of the sum is: %.4f\n', densityCs)
fprintf('The density of the product is: %.4f\n', densityCp)



