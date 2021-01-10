function [num,den]=qMinreal(num,den)
    NumK=num(find(num~=0,1));
    DenK=den(find(den~=0,1));
    z=roots(num);
    p=roots(den);
    [mz,nz] = size(z);
    [mp,np] = size(p);

    % Throw away infinities from zeros
    z = z(isfinite(z));
    mz = length(z);
    mp = length(p);
    iz = ones(mz,1);
    % Loop zeros and match poles
    for i=1:mz
        zi = z(i);
            %tol = 10*abs(zi)*sqrt(eps);
            %tol = abs(zi)*10^-5;
            tol = abs(zi)*10^-4;
        kk = find(abs(p-zi) <= tol);
        if all(size(kk))
            p(kk(1)) = [];
            iz(i) = 0;
        end
    end
    % Eliminate matches in zeros:
    z = z(logical(iz));
    num=(NumK/DenK)*(poly(z));
    den=poly(p);
end
