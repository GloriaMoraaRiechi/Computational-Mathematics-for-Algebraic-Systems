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
A2 = randi(100,7,6);  %Generates a matrix of random integers between 1 and 100 inclusive with 7 rows and 6 columns
B2 = sort(A2, 2);        %Sorts the rows in the matrix in ascending order
disp('The matrix after sorting each row is : ');
disp(B2);

b2 = sort(A2, 1);
disp('The matrix after sorting each column is :')
disp(b2);

B2d = sort(A2, 2, "descend")
disp('The matrix after sorting each row in descending order is :');
disp(B2d);

b2d = sort(A2, 2, "descend")
disp('The matrix after sorting each column in descending order is :');
disp(b2d);

%% Given a mxn matrix A3, write a script that returns the matrix B3 made of rows whose means are smaller of the overall matrix mean of A3
%first claculate the overall mean
%Then calculate the mean of all rows
%Find the rows whose mean is less than that of the matrix
%Display the matrix B of rows whose mean is less
clc
A3 = randi(200, 8, 8);
meanOverall = mean(A3, "all"); %Computes the mean of all the elments in the matrix
meanRows = mean(A3, 2);        %Computes the mean of all rows
rowslessOverallMean = meanRows < meanOverall;
B3 = A3(rowslessOverallMean, :); 
disp('The matrix made of rows whose mean is less than the overall mean is :');
disp(B3);


%% Given a mxn matxix A made of integers, write a script that returns the matrix made of rows of an integer matrix A such that every row returned has a multiple of 7 and a multiple of 5
%First iterate each row of the matrix
%Check each row to see if it contains at least one multiple of 5 and one
%multiple of 7
%Store rows that satisfy both conditions in a new matrix
clc
A = randi(100, 8, 9);              %Generates a random 8x9 matrix with elements from 1 to 100
[m, n] = size(A);
result = ones(m, n);              %Realocation of the result matrix to maximum possible size
resultCount = 0;                   %Initialize counter for the number of rows added

for i = 1:m                        %iterates overs a sequence from 1 to 8
    row = A(i, :);                 %Extracts the i-th row from a matrix. : means take all the elements in that row
    mul5 = any(mod(row, 5) == 0);  %Checks if a row contains elements in the matrix which are multiples of five which is then assigned to mul5
    mul7 = any(mod(row, 7) == 0);  %Checks if a row contains elements in the matrix which are multiples of seven which is then assigned to mul7

    if mul5 && mul7
        resultCount = resultCount + 1;  %Increments the counter
        if resultCount < m
            result(resultCount, :);
        else
            error('The resultCount exceeds the preallocated matrix size')
        end
    end
end

if resultCount > 0
    result = result(1:resultCount, :);
else
    result = [];  % If no rows meet the criteria, return an empty matrix
end
 
disp("The matrix that contains rows which have elements that are both multiple of five and seven is :");
disp(mul5and7);











