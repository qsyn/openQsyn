function C = parallel(A,B)
%parallel returns a QCTRL parrallel connection between two qctrl elements
%or qctrl element and numeruc element
%
% Usage:
%
% ab_parallel = qctrl.parallel(A,B)  returns a parallel connection between
% A and B.

if isa(A,'qctrl') && isa(B,'qctrl')
    if A.sampleTime ~= B.sampleTime
        error('sample time must agree')
    end
    [a,b]=tfdata(A);
    [c,d]=tfdata(B);
    C = locParallel(a,b,c,d);
    Ts = A.sampleTime;
elseif isa(A,'qctrl') && isnumeric(B)
    [a,b]=tfdata(A);
    C = locParallel(a,b,B,1);
    Ts = A.sampleTime;
elseif isnumeric(A) && isa(B,'qctrl')
    [c,d]=tfdata(B);
    C = locParallel(A,1,c,d);
    Ts = B.sampleTime;
else
    error('illigal inputs!')
end

C.sampleTime = Ts;


end

function C = locParallel(a,b,c,d)

ad = conv(a,d);
bc = conv(b,c);
n = max([length(ad) length(bc)]);
if length(ad)<n
    ad = [zeros(1,n-length(ad)) ad];
end
if length(bc)<n
    bc = [zeros(1,n-length(bc)) bc];
end
num = ad+bc;
den = conv(b,d);

Z = roots(num);
P = roots(den);

C0 = qctrl(Z,P,1);

if polyval(num,0)==0 || polyval(den,0)==0
    K1 = 10^(imag(nicresp(C0,1))/20);
    C = C0*K1*polyval(num,1)/polyval(den,1);
else
    C = C0*polyval(num,0)/polyval(den,0);
end


end

