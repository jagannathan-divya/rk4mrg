% generate L matrix in the evolution equation for the coefficient vector

function [L] = Lmat(N)
    n = 0:2:2*(N-1);
    leadDiag = 0.5*(2*n+1);
    upleadDiag = 0.5*sqrt((n(1:end-1)+1).*(n(1:end-1)+2));
    downleadDiag = 0.5*sqrt((n(2:end)).*(n(2:end)-1));
    A0 = diag(leadDiag);
    A1 = diag(upleadDiag,1);
    Am1 = diag(downleadDiag,-1);
    L = A0 + A1 + Am1;
end