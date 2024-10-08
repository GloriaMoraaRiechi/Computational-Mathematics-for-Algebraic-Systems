function b = readMatrixCCW(a)
% SYNTAX: b = readMatrixCCW(a)
% PURPOSE: to read a matrix spiraling clockwise starting
%           with the northwest corner and returning the elements
%           in a row
% INPUT: a = an nxm matrix
% OUTPUT: b = a 1 x nm row

% I need to keep track of row or column, left or right if row,
% up or down if column, beginning and end

% Assume that r1 <= r2  and c1<=c2

[m,n] = size(a);

% a typical step:
% if row r going right, pick the elements, switch to column mode, and the
% next row will be = end +1-r

% Need 4 counters r1, r2, c1, c2
r1 = 1;  r2 = m ; c1 = 1 ; c2 = n;
b= [] ;

rowcol = "row" ;  direction = "right" ;  % clockwise from top left corner

% Main loop
while r1 <= r2 && c1 <= c2
    e = pickelements(a, r1, r2, c1, c2, rowcol, direction) ;
    b = [b, e] ;
    
    
    %update row/column and direction
    if rowcol == "row" & direction == "right"
        rowcol  = "col" ; direction = "down" ;
        r1 = r1+1 ;
    elseif rowcol == "col" && direction == "down"
        rowcol = "row"; 
        direction = "left"; 
        c2 = c2 - 1;  % move the right boundary left
    elseif rowcol == "row" && direction == "left"
        rowcol = "col"; 
        direction = "up"; 
        r2 = r2 - 1;  % move the bottom boundary up
    elseif rowcol == "col" && direction == "up"
        rowcol = "row"; 
        direction = "right"; 
        c1 = c1 + 1;  % move the left boundary right        
    end
end % end of main loop

end

function e = pickelements(a, r1, r2, c1, c2, rowcol, direction)
% Pick elements depending on whether we're moving along a row or column

% if row r and going right from c1 to c2 need  a(r, c1:c2)
if rowcol == "row" & direction == "right"
    e = a(r1, c1:c2 ) ;
elseif rowcol == "col" & direction == "down"
    e = a(r1:r2, c2) ;
    e = e(:)' ;   %flattern column into a row
elseif rowcol == "row" && direction == "left"
    e = a(r2, c2:-1:c1);
elseif rowcol == "col" && direction == "up"
    e = a(r2:-1:r1, c1); 
    e = e(:)';  % flatten column into a row

end

end

function b = readMatrixCCW(a)


% Run learner solution.
a = reshape(1:9, 3, 3) 
b = readMatrixCCW(a) 

% Run reference solution.
bRef = reference.readMatrixCCW(a); 

% Compare.
assessVariableEqual('b', bRef);