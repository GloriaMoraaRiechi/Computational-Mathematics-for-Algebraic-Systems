clc
T=randi(4,1,3)>2;

v=[];
rng=25;
for i=1:100;
    a = toeplitz([T, zeros(1,7)], zeros(1,10));
    v = [v, sum(a(:))];
end
densityRatio = nnz(v)/numel(v)


%% 
clc
% Create the Toeplitz matrix
rng(25); % Set seed for reproducibility
a = toeplitz([randi(4, 1, 3) > 2, zeros(1, 7)], zeros(1, 10));

% Count the number of non-zero elements
num_non_zero = nnz(a); % Count non-zero elements

% Calculate the total number of elements
total_elements = numel(a); % Total number of elements

% Calculate the density ratio
density_ratio = num_non_zero / total_elements; 

% Display the density ratio
disp(['Density Ratio: ', num2str(density_ratio)]);

