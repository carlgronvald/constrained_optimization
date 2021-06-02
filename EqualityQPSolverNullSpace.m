function [x,lambda,time_N] = EqualityQPSolverNullSpace(H,g,A,b)
    [n,m] = size(A);
    start = cputime;
    [Q,R] = qr(A,'vector');
    time_N = cputime-start;
    Qrange = Q(:,1:m);
    Qnull = Q(:,m+1:n);
    R = R(1:m,1:m);
    Y = (R'\b);
    Qnt = Qnull';
    lpre = Qnt*H*Qnull;
    L = chol(lpre);
    mu=L'\(-Qnt*(H*Qrange*Y+g));
    Z=L\mu;
    x = Qrange*Y+Qnull*Z;
    lambda = R\Qrange'*(g+H*x);

    
    