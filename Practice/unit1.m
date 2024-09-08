s1=int2str(3);
s2=num2str(pi);
s=[s1, s2];
num2str(1+2.000)
num2str(1+2.000, '%.3f'); %specifies that the number should be displayed with three digits after the decimal point
clc
% ctr+a shortcut to choose all, ctr-r shortcut to comment at the beginning
% of a line, ctr+t shortcut to remove comment symbol

%% formatSpec
% refres to the formating string used in functions such as fprintf, sprintf
% and num2str to ocntro how the output is dispalyed.
clc
num=123.456789;

str1=num2str(num, '%.3f'); %fixed point with three decimal places
str2=num2str(num, '%d');    %integer format
str3=num2str(num, '%e');   %scientific notation

%% CONDITIONAL AND LOGICAL EXPRESSIONS
clc
if 1==2      %checks if 1 is equal to 2. since it is false, matlab will skip to the next block
    x=3
elseif 2==5  %checks whether 2 is equal to 5, it will skip to the next block since it is false
    x=4
else         %since neither of the above is true, matlab now executes the code in else block
    x=7
end
%% WHITESPACES INSIDE MATLAB ARRAYS
%Refer to any character or sequence of characters that are used to create
%space in a text but are not visible they include;
% spaces
% tabs(\t)
% newlines(\n)
% carriage returns(\r)
% form feeds(\f)
clc
x = input('Enter a number: ');  % Single space around the equal sign for clarity

if x > 5
    fprintf('%d is greater than 5\n', x)  % Indentation improves readability
end
%% ERRORS
x = input('Enter a number: ');

% Intentional runtime error
result = 10 / x;  % If x is 0, this will cause a division by zero error.
disp(['Result is ', num2str(result)]);

%HANDLING RUNTIME ERRORS

try
    result = 0 / x;
catch ME
    disp('An error occurred:');
    disp(ME.message);
end
%% LOGICAL ERRORS

clc
numbers = [1, 2, 3, 4, 5];
total = 0; %initialized to zero and will be used to accumulate the sum of the numbers in the numbers vector

for i = 1:length(numbers)  %this loop iterates over each index of the numbrs vector
    total = total + numbers(i); %in each iteration, thee current nuumber(i) is added to total which will hold the sum of all numbers after the loop completes
end

average1 = total / (length(numbers) - 1);  % Logical error: should divide by length(numbers)
disp(['The average is ', num2str(average1)]);

average2 = total / length(numbers);
disp(['The average is ', num2str(average2)]); %square brackets are used for concatenation of strings
%% TRY-CATCH ONSTRUCT

clc
try
    if a==3
        a=2;
    end
catch
    disp('Unrecognized function or variable a');
    a=9;
end
b=5;

%% TRY-CATCH CONSTRUCT
clc
a1=[2,3];
d=4;
g=1;
    
try                %put in code in a block thet might throw an error
    c=[1, 2] + a1; %Will generate an error if a1 is not 1x2
    d=d+1;         %Will generate an error if d does not exist
    if g==1        %A sppecial error not buil-in to matlab
        errorId = 'Matlab:gIs1'; %MATLAB: is required
        errorMessage = 'g is 1 and it should not be';
        me1 = MException(errorId, errorMessage); %Creates a me1 object
        throw (me1);
    end
    h = 100;       %Will be skipped if any error occurred above
catch me1          %Error in the try-clause will create the me1 object
    disp(me1)
end
%% SYMMETRY OF THE FOR LOOP
%for loop is consider optimistic because the number of iterations is known
%and the loop may break even earlier
clc
xmax=30; n=6; nn=1:n;
x1=1; y1=[];
for i = nn   %i is a member of nn
    disp(i)  %displays the current value of i
    y1 = [y1, x1];  %Appends the current value of x1 t y1
    x1 = 2*x1 + 3;  %Update the valu of x1
    if (x1 > xmax)  %Checks if the update value exceeds xmax
        break       %Terminate the for loop if the updated value of x1 exceeds xmax
    end
end
%disp(sprintf('x is %g \n', y1))
disp(fprintf('X1 is %g \n', y1))  %fprintf returns the number of characters written

% %g automatically selects the most compact form between fixed point
% notation and exponential notation of a given number

%% SYMMETRY FOR WHILE LOOP
%while loop is considered pessimistic because the number of itertaions is
%not known but it is limited to a given number
clc
x2=1; y2=[]; j=1; xmax1=30; n=6; nn=1:n;
while ~(x2 > xmax1)
    disp(j)
    y2 = [y2, x2]
    x2 = 2*x2 + 3;
    j = j + 1;
    if ~ismember(j, nn)
        break
    end
    fprintf('x is %g \n', y2)

end



