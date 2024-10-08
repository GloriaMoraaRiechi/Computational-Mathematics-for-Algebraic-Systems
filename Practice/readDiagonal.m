function diagonalElements = readDiagonal(A)
    [m, n] = size(A);  % Get the size of the matrix
    diagonalElements = [];  % Initialize an empty array to store diagonal elements

    % Traverse each diagonal starting from the first row
    for col = 1:n
        r = 1;
        c = col;
        while r <= m && c > 0
            diagonalElements(end+1) = A(r, c);  % Append elements from the diagonal
            r = r + 1;  % Move down the row
            c = c - 1;  % Move left in the column
        end
    end

    % Traverse each diagonal starting from the second row
    for row = 2:m
        r = row;
        c = n;
        while r <= m && c > 0
            diagonalElements(end+1) = A(r, c);  % Append elements from the diagonal
            r = r + 1;  % Move down the row
            c = c - 1;  % Move left in the column
        end
    end
end
