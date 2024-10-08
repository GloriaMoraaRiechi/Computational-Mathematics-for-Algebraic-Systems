v = []; 
density = []; 
rng(25); % Set the random number generator seed for reproducibility

for i = 1:1000
    a = toeplitz([randi(4, 1, 3) > 2, zeros(1, 7)], zeros(1, 10));
    non_zero_count = sum(a(:) ~= 0); % Count non-zero elements
    v = [v, sum(a(:))]; % Store the sum
    density = [density, non_zero_count / 100]; % Calculate and store density
end

min_density = min(density); % Minimum density across all iterations
disp(min_density); % Display minimum density
