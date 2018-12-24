

k=qpar('k',2,2,5,8);
a=qpar('a',3,1,3,8);
z=qpar('z',0.6,0.3,0.6,8);
wn=qpar('wn',4,4,8,8);

num = [k*wn*wn k*a*wn];
den = [1 2*z*wn wn*wn];
P = qplant(num,den)

w = [0.2 0.5 1 2 5 10 20 50];
cgrid(P,w)
