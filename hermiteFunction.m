% returns a column vector of first n-even modes hermite functions at x
% CAUTION: this is the complete basis iff x=0. For any x \neq 0, you are
% not including the odd modes!!!

function hvec = hermiteFunction(N, x)
    hvec = linspace(0,0,N); % row vector
    hvec(1) = pi^(-.25)*exp(-0.5*x^2);
    iter = 2;
    while iter<=N
            n = 2*(iter-1);
            a = sqrt(2/n); b = sqrt((n-1)/n);
            hvec(iter) = 0 - b*hvec(iter-1);   
            iter = iter + 1;
    end
    hvec = hvec'; %column vector
end

   
    
