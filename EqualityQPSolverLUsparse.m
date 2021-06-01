function [x, lambda] = EqualityQPSolverLUsparse(H,g,A,b)
    KKT = get_KKT_sparse(H,g,A,b);
    [L,U,P] = lu(KKT);
    solution = U \ (L \ (P * -[g;b]));
    x = solution(1:size(H,1));
    lambda = solution(size(H,1)+1:size(H,1)+size(b,1));
end
