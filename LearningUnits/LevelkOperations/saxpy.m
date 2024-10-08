%stands for Single-Precision A-X Plus Y. it is an a Level-1 basic liinear
%algebra subprogram operations. It modifies each element of the vector y by
%adding the corresponding element of the vector x multiplied by the scalar
%alpha. Total flops is 2n for a vector of length n because for each of the
%n elements we perform 1 multiplication and 1 addition

clc
% y = alpha * x + y;
scalar = 2;
x = [1, 2, 3];
y = [4, 5, 6];

%Update each element of y based on the rule provided
for i = 1:length(x)
    y(i)= scalar * x(i) + (4 + (1-i));
end
disp(y)

