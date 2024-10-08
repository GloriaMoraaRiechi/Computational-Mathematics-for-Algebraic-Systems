%% order known
clc

A = input('Enter the 4 by 4 matrix: ')

%Compute the powers
A2 = A^2;
disp('A^2 = ');
disp(A2);
A3 = A^3;
disp('A^2 = ');
disp(A2);
A4 = A^4;
disp('A^2 = ');
disp(A2);

%Check if A^4 is the zero matrix
if isequal(A4, zeros(4))
    disp('The matrix is nilpotent with an index of 4');
else
    disp('The matrix is not nilpotent');
end

%% order not known
clc
A = input('Enter the matrix: ');

%Initialize variables
n = size(A, 1);     %size of the matrix (number of rows)
k = 1;                 %power of the matrix
Ak = A;               %first power of the matrix

%Loop to find the smallest k such that A^k=0
while ~isequal(Ak, zeros(n))
    k = k + 1;
    Ak = A^k;

    %stop if k exceeds n (index of nilpotency can't exceed matrix size)
    if k > n
        break;
    end
end

%Check if nilpotency was found
if isequal(Ak, zeros(n))
    fprintf('The matrix is nilpotent with an index of %d\n', k)
else
    fprintf('The matrix is not nilpotent, checked upto k = %d\n', k)
end
