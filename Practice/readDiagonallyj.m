function diagonals = readDiagonally(A)
    [rows, cols] = size(A);
    diagonals = {}; % Initialize cell array to store diagonals
    
    % Extract diagonals starting from each element in the first row
    for start_col = 1:cols
        r = 1; % Start from the first row
        c = start_col; % Start from the current column
        diagonal = [];
        
        % Collect elements for the current diagonal
        while r <= rows && c >= 1
            diagonal = [diagonal, A(r, c)];
            r = r + 1;
            c = c - 1;
        end
        
        % Store the diagonal in the cell array
        diagonals{end+1} = diagonal;
    end
    
    % Extract diagonals starting from each element in the first column (excluding the top-right corner)
    for start_row = 2:rows
        r = start_row; % Start from the current row
        c = cols; % Start from the last column
        diagonal = [];
        
        % Collect elements for the current diagonal
        while r <= rows && c >= 1
            diagonal = [diagonal, A(r, c)];
            r = r + 1;
            c = c - 1;
        end
        
        % Store the diagonal in the cell array
        diagonals{end+1} = diagonal;
    end
end
