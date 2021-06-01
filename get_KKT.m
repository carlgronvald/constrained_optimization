function [KKT] = get_KKT(H,g,A,b)
    KKT = [H -A;-A', zeros(size(A,2), size(A,2))];