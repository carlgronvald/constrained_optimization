function [x, lambda] = EqualityQPSolverLUdense(H, g, A, b)
    KKT = get_KKT(H,g,A,b);
    [L,U] = lu(KKT);
    solution = U \ (L \ ( -[g;b]));
    x = solution(1:size(H,1));
    lambda = solution(size(H,1)+1:size(H,1)+size(b,1));
end

