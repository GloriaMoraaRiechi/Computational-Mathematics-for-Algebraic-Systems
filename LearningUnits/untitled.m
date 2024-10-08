% Coefficient matrix
A = [2, 1, -1; -3, -1, 2; -2, 1, 2];
% Right-hand side
b = [8; -11; -3];

% Augmented matrix
Ab = [A b];

n = size(Ab, 1); % Number of rows

% Forward elimination (creating the upper triangular form)
for i = 1:n-1
    for j = i+1:n
        factor = Ab(j,i) / Ab(i,i);
        Ab(j, :) = Ab(j, :) - factor * Ab(i, :);
    end
end

disp('Augmented matrix after forward elimination:');
disp(Ab);

% Back substitution
x = zeros(n,1); % Solution vector initialized
x(n) = Ab(n,end) / Ab(n,n);
for i = n-1:-1:1
    x(i) = (Ab(i,end) - Ab(i,i+1:n) * x(i+1:n)) / Ab(i,i);
end

% Display solution matrix
disp('The solution matrix [x1; x2; x3] is:');
disp(x);
%% 

clc
A = [2, 1, -1; -3, -1, 2; -2, 1, 2];
b = [8; -11; -3];
Ab = [A b]; % Augmented matrix
R = rref(Ab); % Reduced Row Echelon Form

% Display the RREF of the augmented matrix
disp('The Reduced Row Echelon Form (RREF) of the augmented matrix is:');
disp(R);

% The last column of the RREF matrix contains the solutions
disp('The solution matrix [x1; x2; x3] is:');
disp(R(:, end));

