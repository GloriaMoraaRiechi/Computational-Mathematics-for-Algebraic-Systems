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
for k = 1 : ...
    e = pickelements(a, r1, r2, c1, c2, rowcol, direction);
    b = [b, e] ;
...
    if rowcol == "row" & direction == "right"
        rowcol  = "col" ; direction = "down" ;
        r1 = r1+1 ;
    elseif 
...
    end
end % end of main loop

end

function e = pickelements(a, r1, r2, c1, c2, rowcol, direction)

% if row r and going right from c1 to c2 need  a(r, c1:c2)
if rowcol == "row" & direction == "right"
    e = a(r1, c1:c2 ) ;
elseif rowcol == "col" & direction == "down"
    e = a(r1:r2, c2) ; e = e(:)' ;
...
end

end

