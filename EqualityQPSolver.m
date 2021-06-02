function [x, lambda,time] = EqualityQPSolver(H, g, A, b, solver)
%EQUALITYQPSOLVER Solves the EQP given by H, g, A, b
%   minimizes 1/2 x' H x + g'x st. A'x=b using the given solver
%
%   Allowed solvers are LDLdense, LDLsparse, LUdense, LUsparse, rangespace
%   & nullspace
time = 0;
if solver == "LDLdense"
    [x, lambda] = EqualityQPSolverLDLdense(H,g,A,b);
elseif solver == "LDLsparse"
    [x, lambda] = EqualityQPSolverLDLsparse(H,g,A,b);
elseif solver == "LUdense"
    [x, lambda] = EqualityQPSolverLUdense(H,g,A,b);
elseif solver == "LUsparse"
    [x, lambda] = EqualityQPSolverLUsparse(H,g,A,b);
elseif solver == "rangespace"
    [x, lambda,time] = EqualityQPSolverRangeSpace(H,g,A,b);
elseif solver == "nullspace"
    [x, lambda,time] = EqualityQPSolverNullSpace(H,g,A,b);
else
    error("solver " + solver + "does not exist; possible values are LDLdense, LDLsparse, LUdense, LUsparse, rangespace, and nullspace")
end

