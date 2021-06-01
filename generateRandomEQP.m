function [H,g,A,b,x,lambda] = generateRandomEQP(n,m)

H = rand(n); %Random n by n matrix
H = 0.5*(H+H') + n*eye(n); %Makes H symmetric and pos def


A = 10*rand(n); 
A = 0.5*(A+A')+n*eye(n);
A = A(:,1:m); 



x = rand(n,1);
lambda = rand(m,1);

KKT = [H -A;-A' zeros(m)];
sol = [x; lambda];
rhs = KKT*sol;

g = -rhs(1:n);
b = -rhs(n+1:n+m);
end







