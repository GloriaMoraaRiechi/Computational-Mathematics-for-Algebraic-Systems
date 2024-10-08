function b = readMatrixCCW(a)
    % SYNTAX: b = readMatrixCCW(a)
    % PURPOSE: to read a matrix spiraling clockwise starting
    %           with the northwest corner and returning the elements
    %           in a row
    % INPUT: a = an nxm matrix
    % OUTPUT: b = a 1 x nm row
    
    [m, n] = size(a);  % Get matrix size
    
    % Initialize row and column boundaries
    r1 = 1; r2 = m;
    c1 = 1; c2 = n;
    
    % Initialize output array and loop control variables
    b = [];
    rowcol = "row";   % row or column
    direction = "right";  % Movement direction
    
    % Main loop continues until all rows/columns are exhausted
    while r1 <= r2 && c1 <= c2
        % Pick elements based on the current direction and orientation
        e = pickelements(a, r1, r2, c1, c2, rowcol, direction);
        b = [b, e];  % Append the picked elements to the result
        
        % Switch direction and update boundaries
        if rowcol == "row" && direction == "right"
            rowcol = "col"; direction = "down";
            r1 = r1 + 1;  % Move top boundary down
        elseif rowcol == "col" && direction == "down"
            rowcol = "row"; direction = "left";
            c2 = c2 - 1;  % Move right boundary left
        elseif rowcol == "row" && direction == "left"
            rowcol = "col"; direction = "up";
            r2 = r2 - 1;  % Move bottom boundary up
        elseif rowcol == "col" && direction == "up"
            rowcol = "row"; direction = "right";
            c1 = c1 + 1;  % Move left boundary right
        end
    end
end  % This 'end' is added to close the main function

function e = pickelements(a, r1, r2, c1, c2, rowcol, direction)
    if rowcol == "row" && direction == "right"
        e = a(r1, c1:c2);  % Pick the top row, moving right
    elseif rowcol == "col" && direction == "down"
        e = a(r1:r2, c2); e = e(:)';  % Pick the right column, moving down
    elseif rowcol == "row" && direction == "left"
        e = a(r2, c2:-1:c1);  % Pick the bottom row, moving left
    elseif rowcol == "col" && direction == "up"
        e = a(r2:-1:r1, c1); e = e(:)';  % Pick the left column, moving up
    end
end


% Run learner solution.
a = reshape(1:9, 3, 3) 
b = readMatrixCCW(a) 

% Run reference solution.
bRef = reference.readMatrixCCW(a); 

% Compare.
assessVariableEqual('b', bRef);