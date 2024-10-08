A = [1 4 7 10; 2 5 8 11; 3 6 9 12]
B = reshape(A, 2,6);
C = reshape(A, 2, 2, 3);
B1 = flipud(A);
B2 = fliplr(A)

%% 

%Shifting and rotating
%This is achieved using the circshift funtion
As = [1 2 3 4; 5 6 7 8; 9 10 11 12];
Bs = circshift(As, [0 2]);; %Shifts the columns to the right by 2 and the rows by 0
Cs = circshift(As, [-1, 0]); %Shifts the rows up by 1 and keeps the columns in place

%Rotating a matrix
%The rot90 function is used to rotate a matrix counterclockwise by 90
%degrees

Br = rot90(As);
Br3 = rot90(Br, 3);

D = magic(6);
Ds = sortrows(D);  %Sorts in ascending oredr according to the elements in the firts column
%% Creating Multidimensional Arrays
clc
A = [1 2 3; 4 5 6; 7 8 9];
A(:,:,2) = [10 11 12; 13 14 15; 16 17 18];
%The cat function can be a useful tool for building multidimensional arrays. For example, create a new
%3-D array B by concatenating A with a third page. The first argument indicates which dimension to
%concatenate along.
B = cat(3,A,[3 2 1; 0 9 8; 5 3 7]);

%Another way to quickly expand a multidimensional array is by assigning a single element to an entire
%page. For example, add a fourth page to B that contains all zeros.

B(:,:,4) = 0

elB = B(3,3,3); %Accessing the element in the 3rd page, 3rd row and 3rd column

C = B(:, [1 3], :); %Accessing the first and last columns of each page

%find the second and third rows of each page, use the colon operator to create your index vector.
D = B([2 3], :, :);

%% Manipulating Arrays
%This can be achieved using reshape, permute and squeeze functions which
%are useful for rearranging elements
clc
A = [1 2 3 4 5; 9 0 6 3 7; 8 1 5 0 2];
A(:,:,2) = [9 7 8 5 2; 3 5 8 5 1; 6 9 4 3 3]

%The reshape function operates clockwise creating a new matrix by taking
%consecutive elements down each column of A starting with the first spage
%and then moving to the second page
B = reshape(A, [6 5]);

%permutations are used to rearrange the order of the dimensions of an array
P1 = permute(A,[2 1 3]);


%% Squeeze function

%When working with multidimensional arrays, you might encounter one that has an unnecessary
%dimension of length 1. The squeeze function performs another type of manipulation that eliminates
%dimensions of length 1. For example, use the repmat function to create a 2-by-3-by-1-by-4 array
%whose elements are each 5, and whose third dimension has length 1.
clc
A = repmat(5,[2 3 1 4]);
size(A);
numdimsA = ndims(A);

%% 
%Use the squeeze function to remove the third dimension, resulting in a 3-D array

B = squeeze(A);
size(B);
numdimsB = ndims(B);