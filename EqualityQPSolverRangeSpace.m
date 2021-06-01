function [x, lambda] = EqualityQPSolverRangeSpace(H,g,A,b)
    M=chol(H); 
    mu=M'\g;
    Hg=M\mu;
    mu=M'\A;
    HA=M\mu;
    lambda = (A'*HA)\(b+A'*Hg);
    x = HA*lambda-Hg;