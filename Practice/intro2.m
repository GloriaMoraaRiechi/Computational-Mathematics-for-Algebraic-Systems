clc
a = 6*(10/13) + 18/(5*7) + 5*(9^2);
6*(35^(1/4)) + 14^0.35;

%% 
%finding the radius of a circular cylinder whose volume is 20% greater of
%another cylinder of the same height
clc
r = 8;
h = 15;
V = pi*r^2*h;
V = V + 0.2*V;
r = sqrt(V/(pi*h));
r
%% 
clc
x=-5+9i;
y=6-2i;
x+y;
x*y;
x/y;

u=0:0.1:1.0;
w=5*sin(u);
m=length(w);
%% polynomial roots
%polynomials can be described in matlab as as arrat whose elements are the
%polynomial coefficients starting with the coefficient of the highest power
%of x> the roots of a polynomial f(x) are the values such that f(x)=0 and
%can be found with the roots(a) function

%find the roots of x^3-7x^2+40x-34=0
clc
p1=[1 -7 40 -34];
roots(p1); 
roots([1 -7 40 -34]);

%find the roots of 290-11x+6x^2+x^3
%first we need to rearrange the polynomial to x^3+6x^2-11x+290
p2=[1 6 -11 290];
roots(p2)
%% 
clc
m1=cos(0):0.02:log10(100);
length(m1);
m1(25);

sin(90);
sind(90);
%% PLOTTING WITH MATLAB
clc
x1=0:.01:7;
y1=3*cos(2*x1);

plot(x1, y1), xlabel('x1'), ylabel('y1')
title('Graph of y=3cos(2x)')
%% OVERLAY PLOTS
clc
x2=0:0.01:5;
y2=2*sqrt(x2);
z2=4*sin(3*x2);
plot(x2, y2, x2, z2, '--');
xlabel('x2'), gtext('y2'), gtext('z2');
%% CREATING LINE PLOT FROM MATRIX
clc
y3=magic(4);
plot(y3);
%% SPECIFY LINE STYLE
x4=0:pi/100:2*pi;
y4=sin(x4);
y5=sin(x4-0.25);
y6=sin(x4-0.5);

figure
plot(x4, y4, x4, y5, '--', x4, y6, ':');
%% SPECIFY LINE STYLE, COLOR AND MARKER
clc
x5=0:pi/10:2*pi;
y7=sin(x5);
y8=sin(x5-0.25);
y9=sin(x5-0.5);


%plot(x5, y7, 'g')  %Green normal line with no markers
%plot(x5, y8, 'b--o') %a blue dashed line with circle markers
%plot(x5, y9, 'c*')  %cyyan stars 

plot(x5,y7,'g',x5,y8,'b--o', x5,y9,'c*')
xlabel('x5'), gtext('y7'), gtext('y8'), gtext('y9')

%% DISPLAY MARKERS AT SPECIFIC DATA POINTS
%create a line plot and display markers at every fifth data point by
%specifying a marker symbol and setting the MarkerIndices property as a
%name-value pair
clc
x6=linspace(0,10,200); %creates a vector with 100 points linearly spaced from 0 to 6
y10=sin(x6);
plot(x6,y10, '-o', 'MarkerIndices',1:55:length(y10)) 
%-o specifies that the plot should use a line with circle markers
%MarkerIndices 1:5:length(y10) specifies that the markesr should be placed
%every 55 indices along the length of y10 starting from the fisrt index.
%This places markers at approximately every 55th point in the plot
%% SPECIFY LINE WIDTH, MARKER SIZE AND MARKER COLOR
clc
x11=-pi:pi/10:pi;
y11=tan(sin(x11))-sin(tan(x11));
plot(x11,y11,'--gs', ...   %dashed green line with square markers
    'LineWidth',2, ...    %set the width of the line to 2 points
    'MarkerSize', 10, ...  %sets the size of the marker to 10 points
    'MarkerEdgeColor','b', ...  %sets the colors of the edge of the markers to blue
    'MarkerFaceColor',[0.5,0.5,0.5])%sets the color of the inside of the markers to gray
%% plot the function s=2sin(3t+2)+sqrt(5t+1) over the interval 0 to 5.Put a title on the plot and label the axis
%s represents speed in feet per second
%t represents time in seconds


clc
t=0:5;
s=2*sin((3*t)+2) + sqrt((5*t)+1);
plot(t, s, ':m', 'LineWidth',2)
xlabel('Time in seconds'), ylabel('Speed in feet per second')
title('Graph of s=2sin(3t+2)+sqrt(5t+1)')
%% plot the functions y12=4*sqrt(6*x12+1) and z12=5*exp(0.3*x12)-2*x12 over the inter 0 to 1.5

clc
x12=0:.09:1.5;
y12=4*sqrt(6*x12+1);
z12=5*exp(0.3*x12)-2*x12;
plot(x12, y12, x12, z12)
grid on
xlabel('Distance in meters'), ylabel('Force in Newtons')
gtext('y12'), gtext('z12')

%% SUBPLOTS

