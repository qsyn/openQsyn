%% Example:Plant with unstructured uncertainty
%
%%
% The plant is given as 
%
% $$ P(s) = \frac{s+a}{1 + s \zeta / \omega_n + s^2 / \omega_n^2} $$
%  
% with uncertain paraetmers given as 
% 
% $$ k \in [2,5],~ a \in [1,3],~ \zeta \in [0.1,0.6],~ \omega_n \in [4,8]$$
% 
% and with unstructured uncertainty given as 
w = [0.1 0.2 0.5 1 2 5 10 20 50 100];
m = [0 0.3 0.3 0.3 0.3 0.35 0.35 0.35 0.5 0.5];
semilogx(w,m); xlabel('m(j\omega)'); ylabel('freqeuncy [rad/s]')
%%
% Define uncertain parameters:
k=qpar('k',2,2,5,8);
a=qpar('a',3,1,3,8);
z=qpar('z',0.6,0.1,0.6,8);
wn=qpar('wn',4,4,8,8);

%%
% Construct the numerator and denomerator and plant
num = [k*wn*wn k*wn*wn*a];
den = [1 2*z*wn wn*wn];
P = qplant(num,den)

%%
% compute the nominal and template by e.g. recurcive grid ans show the
% results
P.ctpl('recgrid',w);
P.cnom(logspace(-2,3,200));  
P.showtpl
%%
% Add the unstructored uncertainty unsing the command aunstruc - a shortcut
% for add unstructored uncertainty
P.aunstruc(w,m)

%%
% Now compute the templates again 
P.ctpl('recgrid',w);
P.showtpl
