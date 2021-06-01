function [x, lambda] = EqualityQPSolverLDLsparse(H, g, A, b)
    KKT = get_KKT_sparse(H,g,A,b);
    [L,D, P] = ldl(KKT);
    solution = P*(L' \ (D \ (L \ (P'*-[g;b]) )));
    x = solution(1:size(H,1));
    lambda = solution(size(H,1)+1:size(H,1)+size(b,1));
