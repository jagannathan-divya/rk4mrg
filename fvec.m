function f = fvec(N,gmma)
    n = 0:2:2*(N-1);

    % The second derivative matrix
    leadDiag = -0.5*(2*n+1);
    upleadDiag = 0.5*sqrt((n(1:end-1)+1).*(n(1:end-1)+2));
    downleadDiag = 0.5*sqrt((n(2:end)).*(n(2:end)-1));
    D0 = diag(leadDiag);
    D1 = diag(upleadDiag,1);
    Dm1 = diag(downleadDiag,-1);
    Dsq = D0 + D1+ Dm1;

    % Construction of f vector
    f = zeros(N,1);
    iter = 1;
    while iter<=N
        % indicator vector
        ivec = zeros(N,1); ivec(iter) = sqrt(2*pi)*((-1)^(iter-1))*(gmma/pi);
        % di vector
        di = (gmma^2*eye(N)-transpose(Dsq))\ivec;
        f(iter) = dot(di,hermiteFunction(N,0));
        iter = iter + 1;
    end
end