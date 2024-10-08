%% 

function seq = readDiagonals(A)

    [m, n] = size(A);  % Get the dimensions of the matrix
    seq = [];          % Initialize the output sequence

    % Traverse diagonals starting from top-right to bottom-left
    for startCol = n:-1:1
        row = 1;
        col = startCol;
        while row <= m && col >= 1
            seq = [seq, A(row, col)];
            row = row + 1;
            col = col - 1;
        end
    end
    
    % Traverse diagonals starting from the columns of the first row (excluding the first column)
    for startRow = 2:m
        row = startRow;
        col = n;
        while row <= m && col >= 1
            seq = [seq, A(row, col)];
            row = row + 1;
            col = col - 1;
        end
    end


%% 
clc
A = randi(10, 11)
Ba = A(1:5, 5:10);

Bb = A(6:10, 6:10);

 Bc = [] ; for i=6:10, Bc = [Bc; A(i, 7:11) ] ; end


Bc

% Bd = [] ; for j=7:11, Bd = [Bd ; A(6:10, j) ] ; end
% Bd


Be = [] ; for j=7:11, Be = [Be, A(6:10, j) ] ; end

Be

Bb

e1 = length(A)*length(A');
e2 = size(A,1)*size(A,2);

s = sum(sum(A<0))
e3 = length(A(:));

[A(:,2) , A(:,3)]';
















