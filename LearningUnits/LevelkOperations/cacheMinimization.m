%Suppose your cache memory has only 5 places that can be read or written in one shot. Suppose you pay $10 every time the cache is activated. what is the cheapest way to compute Ab  when
% A=[1 2 3 4; 5 6 7 8] b=[4;3;2;1]

%STEP1: CACHE LOAD 1
%load the first rwo of A and the entire vector b into the cache

clc
% Initialize matrix A and vector b
A = [1 2 3 4; 5 6 7 8];
b = [4; 3; 2; 1];

Ab = zeros[2, 1];   %holds the output of the matrix-vector multiplication

% Number of cache activations
cache_activation_cost = 0;    %used to track the total cost incurred from activating the cache during calculations

% Compute the dot products row-wise
for i = 1:size(A, 1)
        % Simulate cache load for row i of A and all of b
    cache = [A(i, :), b']; % Load the entire row of A and b
    cache_activation_cost = cache_activation_cost + 10; % Activation cost
    
        % Compute dot product for row i
    result(i) = cache(1) * cache(5) + cache(2) * cache(6) + cache(3) * cache(7) + cache(4) * cache(8);
end

% Display the result
disp('Result of Ab:');
disp(result);

% Display the total cache activation cost
disp('Total cache activation cost:');
disp(cache_activation_cost);





%% 
% Initialize matrix A and vector b
A = [1 2 3 4; 5 6 7 8];
b = [4; 3; 2; 1];

% Initialize the result vector
result = zeros(2, 1);

% Initialize a variable to track the number of cache activations
cache_activation_cost = 0;

% Simulate cache load for the first row of A and the first element of b
cache = [A(1, 1:4), b(1)]; % Load the first row and b(1)
cache_activation_cost = cache_activation_cost + 10;

% Compute partial dot product for the first row
partial_result1 = cache(1) * cache(5); % A(1,1) * b(1)

% Simulate cache load for the remaining elements of b
cache = [A(1, 2:4), b(2:4)]; % Load A(1,2:4) and b(2:4)
cache_activation_cost = cache_activation_cost + 10;

% Complete the dot product for the first row
partial_result1 = partial_result1 + cache(1) * cache(5) + cache(2) * cache(6) + cache(3) * cache(7);
result(1) = partial_result1;

% Simulate cache load for the second row of A and the first element of b
cache = [A(2, 1:4), b(1)];
cache_activation_cost = cache_activation_cost + 10;

% Compute partial dot product for the second row
partial_result2 = cache(1) * cache(5); % A(2,1) * b(1)

% Simulate cache load for the remaining elements of b
cache = [A(2, 2:4), b(2:4)];
cache_activation_cost = cache_activation_cost + 10;

% Complete the dot product for the second row
partial_result2 = partial_result2 + cache(1) * cache(5) + cache(2) * cache(6) + cache(3) * cache(7);
result(2) = partial_result2;

% Display the result
disp('Result of Ab:');
disp(result);

% Display the total cache activation cost
disp('Total cache activation cost:');
disp(cache_activation_cost);

