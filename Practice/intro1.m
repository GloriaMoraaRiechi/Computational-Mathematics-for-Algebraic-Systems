clc
b=[1;2;3;4;5;6;6];

%solve the following simultaneously
%x+y=1 and 2x-5y=9
A=[1 1; 2 -5];
B=[1;9];
%this will be solved using the left division operator
C=A\B;


%POLYNOMIALS
%Polynomials are functions that are built by simply adding together or
%subtracting some power functions. The coefficients of a polynomial are
%entered as a roww vector beginning with the highest power and including
%the ones that are close to zero

%Create a row vector for the following function: y = 2x4 + 3x3 + 5x2 + x + 10

p2=[2 3 5 1 0];

%Create a row vector for the following function: y = 3x4 + 4x2 − 5 

p1=[3 0 4 0 -5];

% The polyval Function
%We can evaluate a polynomial p for a given value of x using the syntax
%polyval(p, x) where p contains the coefficients of polynomial and x is the
%given number

%Evaluate f(x) at 5. f(x)=3x2+2x+1
p3=[3 2 1];
x=5;
polyval(p3,x);

%The roots Function
%In matlab the roots function is used to find the roots easily
%Find the roots for the following: 0.6x2+0.3x-0.9=0
p4=[0.6 0.3 -0.9];
r1=roots(p4);



e1=6*7+4^2-2^4;
e2=(3^2+2^3)/(4^5-5^4) + (64^0.5-5^2)/(4^5+5^6+7^8);
e3=log10(10^2)+10^5;
e4=exp(2)+2^3-log(exp(2));
e5=sin(2*pi)+cos(pi/4);
e6=tan(pi/3)+cos(270*pi/180)+sin(270*pi/180)+cos(pi/3);
p7=[2 4; 1 5];
m7=[1;2];
v7=p7\m7;
y8=polyval([4 0 3  -1 0], 5);

 %% 
 
x=(-pi:.1:pi);
y=sin(x);
plot(x,y);
title('Graph of y=sin(x)');
xlabel('X');
ylabel('sin(x)');
grid on;

%% 

x = 0:.2:12;
Y = [besselj(1,x); besselj(2,x); besselj(3,x)];
plot(x,Y);
legend('First','Second','Third','Location','NorthEastOutside')
b = bar(rand(10,5),'stacked'); 
hold on
ln = plot(1:10,5*rand(10,1),'-o');
hold off
legend([b,ln],'Carrots','Peas','Peppers','Green Beans',...
 'Cucumbers','Eggplant');
 %% SUPERIMPOSED PLOTS

 %This script generates sin(x) and cos(x) plot on the same graph
 %initialization of variables
 x=[-pi:.1:pi]; %creates a row vector from -pi to +pi with increments of 0.1
 y0=sin(x); %calculates sine value for each x
 y1=cos(x); %calculate cosine value for each x

plot(x, y0, x, y1); %plot sin(x) and cos(x) on the same graph
title('Graph of sin(x) and cos(x)'); %title of the graph
xlabel('x'); %label of x axis
ylabel('sin(x), cos(x)'); %label of y axis
legend('sin(x)', 'cos(x)', 'Location', 'best'); %inserts legend in the same order as y0 and y1
grid on; %turns graph grid on
%% MULTIPLE PLOTS IN A FIGURE

%Multiple plots in a figure can be generated with subplot in the command
%time however the built-in plot tools can be used.
close all
%this script will generate sin(x) and cos(x) variables
X1=[-2*pi:.1:2*pi]; %creates a row vector from -2*pi to 2*pi with increments of 0.1
Y1=sin(X1); %calculates sine value for each x
Y2=cos(X1); %claculates cosine value for each x
Y3=Y1+Y2; %calculate sin(x)+cos(x)
Y4=Y1-Y2; %calculates sin(x)-cos(x)
plot(X1, Y1, X1, Y2, X1, Y3, X1, Y4)
grid on