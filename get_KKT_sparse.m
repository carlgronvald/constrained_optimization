function [KKT] = get_KKT_sparse(H,g,A,b)
    KKT = sparse([H -A;-A', zeros(size(A,2), size(A,2))]);