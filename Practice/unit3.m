clc
glyphplot(ones(1, 17)*0.5, 'glyph', 'face', 'standardize', 'off');

%Generate a set of 20 rendom faces using 16 randomly distributed numbers
%(Keep face size, x1 the same
figure
rand('state',43445);

faces1 = 0.1+0.8*rand(20, 17) ;
faces1(:,1) = 0.88*ones(20,1) ;

glyphplot(faces1, 'glyph','face', 'standardize','off')
plot(faces1)
xlabel('features')
ylabel('feature Values')
title('')