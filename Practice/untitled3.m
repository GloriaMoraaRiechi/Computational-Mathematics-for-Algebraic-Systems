A = randi(10,11)
[A(:,2) , A(:,3)]';
[A(3,:) , A(2,:)];
[A(:,2) ; A(:,3)];


e1 = length(A(:))/size(A,1);
e2 = length(A(1,:));
e3 = size(A,2);





