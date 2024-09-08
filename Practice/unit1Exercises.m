H = hilb(6);         % Generate the Hilbert matrix of order 6
H_power = H.^0.1;    % Raise each element of the matrix to the power of 0.1
singular_values = svd(H_power);  % Compute the singular values

% Plotting the singular values
plot(singular_values, 'o-');
title('Singular Values of hilb(6).^0.1');
xlabel('Index');
ylabel('Singular Value');
grid on;

%% Given a mxn matrix A, write a script that sorts every row separately in ascending order

clc
A = randi(100,7,6)  %Generates a matrix of random integers between 1 and 100 inclusive with 7 rows and 6 columns


