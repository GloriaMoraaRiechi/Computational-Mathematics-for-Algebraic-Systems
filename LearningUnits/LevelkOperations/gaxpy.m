%Gaxpy stand for General A.X Plus Y.
%it is a generalization of the saxpy operation used to represent operations
%where a scalar is multiplied by a vector and then added to another vector
%not limited to single precision, can handle various precisions and data
%types. y = alpha*A*x + y where A is a matrix
%can be interpreted as a scaled matrix-vector multiplication followed by a
%vector addition

%Process
%for each element in the vector y, the corresponding row of matrix A is
%multiplied by the vector x, scaled by alpha and then added to the
%corresponding element in the vector y
%FLOPS=m*(2n-1)+m, m*n multiplications and m*(n-1) additions

alpha = 2;
A = [1, 2, 3; 4, 5, 6; 7, 8, 9];
x = [1; 2; 3];
y = [1; 2; 3];
y = alpha * (A * x) + y
