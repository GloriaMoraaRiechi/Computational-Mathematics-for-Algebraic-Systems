% toeplitz( [randi(4,1,3)>2, zeros(1,7)] , zeros(1,10))

v = [ ] ;  rng(25);
for i=1:1000;

            a = toeplitz( [randi(4,1,3)>2, zeros(1,7)] , zeros(1,10));

            v=[v,sum(a(:))];
            density=

end
%% 
clc
v = [];  % Initialize the vector to store densities
rng(25);  % Set random seed for reproducibility
for i = 1:1000
    a = toeplitz([randi(4,1,3)>2, zeros(1,7)], zeros(1,10));
    v = [v, sum(a(:))];  % Record the density (sum of non-zero elements)
end

% Calculate the average density
average_density = mean(v);

% Display the average density
disp(['Average Density: ', num2str(average_density)]);

%% 
clc
v = [];  % Initialize an empty vector to store densities
rng(25);  % Set random seed for reproducibility
for i = 1:1000
    a = toeplitz([randi(4,1,3)>2, zeros(1,7)], zeros(1,10));
    v = [v, sum(a(:))];  % Record the density (sum of non-zero elements)
end

% Display the first element of v
first_element = v(1);
disp(['First element of v: ', num2str(first_element)]);

% Calculate the most frequent matrix density
most_frequent_density = mode(v);
disp(['Most frequent matrix density: ', num2str(most_frequent_density)]);


% Determine the maximum and minimum matrix densities
max_density = max(v)  % Maximum density
min_density = min(v) % Minimum density


%% 
clc
% Set the number of simulations
num_simulations = 1000; 
zero_matrix_count = 0;  % Counter for zero matrices

% Set random seed for reproducibility
rng(25); 

for i = 1:num_simulations
    % Generate the 10x10 Toeplitz matrix
    a = toeplitz([randi(4,1,3) > 2, zeros(1, 7)], zeros(1, 10));
    
    % Check if the matrix is a zero matrix
    if all(a(:) == 0)
        zero_matrix_count = zero_matrix_count + 1;  % Increment the counter
    end
end

% Calculate the probability of generating a zero matrix
probability_zero_matrix = (zero_matrix_count / num_simulations) * 100;  % Convert to percentage

% Display the result
disp(['Probability of having a zero matrix: ', num2str(probability_zero_matrix), '%']);

