function [x, lambda] = EqualityQPSolverLDLdense(H,g,A,b)
    KKT = get_KKT(H,g,A,b);
    [L,D,p] = ldl(KKT,'lower','vector');
    rhs =  -[g;b];
    solution(p) = L' \ (D \ (L \rhs(p)));
    x = solution(1:size(H,1));
    lambda = solution(size(H,1)+1:size(H,1)+size(b,1));


