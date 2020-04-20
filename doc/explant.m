function y = explant(k,a,zet,wn,s)


Ts = 0.02;
z = exp(Ts*s);                              % approximate z

N = length(s);
    
y = zeros(N,1);
for ii=1:N
    num = k(ii).*[1 a(ii)];
    den = conv([1/(wn(ii)*wn(ii)) 2*zet(ii)/wn(ii) 1],[1/10000 2*0.7/100 1]);
    Pz = c2d(tf(num,den),Ts,'zoh');             % convert by ZOH
    [dnum, dden] = tfdata(Pz);                  % cell object are returned
    y(ii) = polyval(dnum{1},z(ii))./polyval(dden{1},z(ii));
end

end

