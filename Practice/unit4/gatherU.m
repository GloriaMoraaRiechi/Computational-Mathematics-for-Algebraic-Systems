function w = GatherU(u)
    % GatherU function to convert a full-storage vector u to sparse format
    % Inputs:
    %   - u: A full-storage vector (can contain zero and nonzero elements)
    % Outputs:
    %   - w: A matrix with two columns: 
    %           - 1st column: Indices of the nonzero elements of u
    %           - 2nd column: Values of the nonzero elements
    
    n = length(u);   % Get the number of elements in u
    nnz = 0;         % Initialize the count for nonzero elements
    
    % Loop through each element of u
    for i = 1:n
        if u(i) ~= 0  % Check if the element is nonzero
            nnz = nnz + 1;       % Increment the nonzero counter
            w(nnz, 1) = i;       % Store the index of the nonzero element
            w(nnz, 2) = u(i);    % Store the value of the nonzero element
        end
    end
    
    % If no nonzero elements were found, return an empty array
    if nnz == 0
        w = [];  % Return an empty array if no nonzero elements were found
    end
end

% Example Usage
u = [0, 5, 0, 0, 3, 0, 7];  % Input vector with zero and nonzero elements

% Call the GatherU function and store the result in w
w = GatherU(u);

% Display the result
disp('Sparse format (index, value):');
disp(w);
