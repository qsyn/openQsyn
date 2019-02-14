function y = explant(k,a,zet,wn,s)

Ts = 0.02;
z = exp(Ts*s);                              % approximate z
num = k*[1 a];
den = conv([1/(wn*wn) 2*zet/wn 1],[1/10000 2*0.7/100 1]);
Pz = c2d(tf(num,den),Ts,'zoh');             % convert by ZOH
[dnum, dden] = tfdata(Pz);                  % cell object are returned
y = polyval(dnum{1},z)./polyval(dden{1},z);

end

