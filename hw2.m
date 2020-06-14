clear 
clc
k = qpar('k',2,2,5,8);
a = qpar('a',3,1,3,8);
z = qpar('z',0.6,0.3,0.8,8);
wn = qpar('wn',4,4,8,8);


num = [k*wn*wn k*a*wn*wn];
den = [1 2*z*wn wn*wn];

P = qplant(num,den);

w_nom = logspace(-2,2,200);
P.cnom(w_nom);


w = [0.4 0.5 1 3 5 7 10 20 50];
P.ctpl('recedge',w);
% P.ctpl('recgrid',w)
spec1 = qspc('odsrs',w,6);
spec2 = qspc.rsrs([1.2 0.2],10,1.5,[],logspace(-1,2),2.85,3.1);

des = qdesign(P,[spec1 spec2]);
des.cbnd('rsrs')
des.cbnd('odsrs')

des.showbnd('odsrs')
des.showbnd('rsrs')
 
h = des.showbnd('odsrs',[],[3 5 7 10 20 50]);
des.showbnd('rsrs',h,[0.4 0.5 1]);

s = zpk(0,[],1);
set(s,'DisplayFormat','Frequency');

% G_Lead = Clead(40,10,1)*Clead(42,30,1);
% G_Lag = Clag(50,1000);

% curr_G= ((s+5)/(s*(s^2 +s +1)));%*G_Lead*G_Lag;
G = qctrl(-1,[0 -1 -2],1);
des.loopnic(G)
ngrid

% spec2.show('freq');
% des.clmag(G,1/(0.5*s+1))
% ylim([-55 10])
%% Filter Design - F
F = 1/(0.5*s+1);
spec2.show('freq');
des.clmag(G,F)
ylim([-55 10])
%% Validation
L = series(P,G); % open loop
S = feedback(L,1); % closed loop from d to y (sensitivity)
T = series(S,series(L,F)); % closed loop from r to y
spec1.show; hold on % show the specs.
pgrid = P.pars.sample(100); % generates 20 random samples
S.bodcases([],w_nom,'showphase',0) % plot magnitude response
xlim([0.01 100])
spec2.show('freq'); hold on % show the specs.
T.bodcases([],w_nom,'showphase',0) % plot magnitude response
axis([0.01 100 -55 10])