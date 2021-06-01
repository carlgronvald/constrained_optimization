function [x,lambda] = EqualityQPSolverNullSpace(H,g,A,b)
    [n,m] = size(A);
    [Q,R] = qr(A);
    Qhat = Q(:,1:m);
    Qn = Q(:,m+1:n);
    R = R(1:m,1:m);
    Y = (R'\b);
    Qnt = Qn';
    Z = (Qnt*H*Qn)\(-Qnt*(H*Qhat*Y+g));
    x = Qhat*Y+Qn*Z;
    lambda = R\Qhat'*(g+H*x);
