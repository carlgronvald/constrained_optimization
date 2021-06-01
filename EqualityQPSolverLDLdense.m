function [x, lambda] = EqualityQPSolverLDLdense(H,g,A,b)
    KKT = get_KKT(H,g,A,b);
    [L,D] = ldl(KKT);
    solution = L' \ (D \ (L \ -[g;b]));
    x = solution(1:size(H,1));
    lambda = solution(size(H,1)+1:size(H,1)+size(b,1));


