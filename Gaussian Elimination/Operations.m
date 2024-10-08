%This is an implementation of Gaussian elimination to solve a system of linear equations A * X = B. 
% Here, A is the coefficient matrix, B is the constant (source) vector, and X is the solution vector.


clc
% A = [10 3 1; 3 10 2; 1 2 20];
% B = [19; 29; 35]

A = input('Enter your coefficient Matrix: ');  %Matrix containing the coefficints of variables in a system
B = input('Enter source vector: ');   %constant terms on the right-hand side of the equations

N = length(B);    %finds the number of equations which is equal to the number of elements in vector B
X = zeros(N,1);   %Initializes the solution vector X with N rows and 1 columns with zeros which will be updated later

aug = [A B];       %Forms the augmented matrix that will allow elimination directly on both the coefficients and the constants

operationCount = 0;

%FORWARD ELIMINATION
%Reduces the system to an upper triiangular matrix where all the elements
%below the diagonal of the matrix are zeros where j is the pivot column
%which controls the column you are currently working on to eliminate the
%elements below the diagonal, N is the number of rows
for j=1:N-1   %the loop starts from the first column and runs upto the second last column
    for i=j+1:N  %for each column j, the loop iterates over the rows i below row j and eliminates the variable in column j from those rows
        m = aug(i,j) / aug(j,j);   %ratio of the current element aug(i,j) and the pivot element aug(j,j) whic is a mutiplier that will be used to eleminate the element aug(i, j)
        aug(i, :) = aug(i, :) - m * aug(j, :);  %updates row i reducing the system to an upper triangular form

        operationCount = operationCount + 1;

        %Display the operation details
        fprintf('Operation %d: Row %d = Row %d - (%.2f) * Row %d\n', operationCount, i, i, m, j);

        % Display the updated augmented matrix after each row operation
        disp('Updated Augmented Matrix:');
        disp(aug);

    end
end
aug;


%BACKWARD SUBSTITUTION
X(N) = aug(N, N+1) / aug(N, N)
for k = N-1: -1:1
    X(k) = (aug(k, N+1) - aug(k, k+1:N) * X(k+1:N)) / aug(k,k);
end
X

% Display the total number of row operations
fprintf('Total number of row operations: %d\n', operationCount);