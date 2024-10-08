clc
%Computes results of non-zero elements
% Define sparse vectors x and y
x = sparse([0; 2; 0; 4]);
y = sparse([1; 0; 3; 0]);

% Define the scalar alpha
alpha = 3;

% Perform sparse SAXPY operation
y = alpha * x + y;

% Display the result
disp(full(y));  % Use full to display the result as a regular vector
