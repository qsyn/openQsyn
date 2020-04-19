%% Example: RFF template computation
% 
% If an uncertain transfer functions is given in Real Factored Form
%
% $$P(s) = k e^{-\tau s} \frac{\displaystyle\prod_{k=1} (1+s/b_{k})}{\displaystyle\prod_{u=1} (1+s/b_{u})}  
%          \frac{\displaystyle\prod_{m=1} (s+b_{m})}{\displaystyle\prod_{v=1} (s+b_{v})} 
%          \frac{\displaystyle\prod_{q=1} (1+2\zeta_q s/\omega_{q} + s^2/\omega_q^2)}{\displaystyle\prod_{w=1} (1+2\zeta_w s/\omega_{w} + s^2/\omega_w^2)}
%          \frac{\displaystyle\prod_{r=1} (s^2+2\zeta_r \omega_{r} s + \omega_r^2)}{\displaystyle\prod_{i=1} (s^2+2\zeta_i \omega_{i} s + \omega_i^2) }
%          \left(1+M(s)\right)
% $$
%
% a method given in Gutman et. al. (1994) can be used to compute accurate templates in a short time.  
% Note that in the above all parameters are uncertain, and $M(s)$ denotes multiplicative unstructured uncertainty
% First and second order factors whose gain equals 1 for s=0 are said to be given in direct current or dc-form, and the
% remaining first and second order factors are said to be in high frequency or hf-form.

%%
% To use the |rff| method, the qplant must be constructed using |qrff|
% elements instead of |qpoly| elements for the numerator and denumerator.

%%
% *Exmaple 1*: 
% 
% The plant is given as 
%
% $$ P(s) = \frac{\displaystyle s+a}{\displaystyle 1 + 2 \zeta s / \omega_n + s^2 / \omega_n^2} $$
%  
% with uncertain paraetmers given as 
% 
% $$ k \in [2,5],~ a \in [1,3],~ \zeta \in [0.1,0.6],~ \omega_n \in [4,8]$$


%%
% We construct the qpar elements as usual
k=qpar('k',2,2,5,8);
a=qpar('a',3,1,3,8);
z=qpar('z',0.6,0.3,0.6,8);
wn=qpar('wn',4,4,8,8);

%%
% Now, we construct num and den as array of |qrff| objects. 
% The numerator has a gain of |k| and a first order hf element; the
% denomerator has a single 2nd order dc element. 
num = [qrff('hf',a) qrff('gain',k)];
den = qrff('dc',wn,z);
P1 = qplant(num,den);

%%
% From this stage we can continue as usual and compute the templates. 
w = [0.2 0.5 1 2 5 10 20 50];
P1.cnom;
P1.ctpl('rff',w,'accuracy',[1 1]);
P1.showtpl

%%
% Note that the fact that the plant is constucted using |qrff| elements
% does not prevent template computations via other methods. For exmaple:

P1.ctpl('random',w,'union',1);
P1.showtpl

%% 
% However, the |rff| method cannot be used on qplants constructed using
% qpoly elements.

%% Example 2: 
% 
% The plant is given as 
%
% $$ P(s) = \frac{\displaystyle k}{\displaystyle s^2} e^{-h s} $$
%  
% with uncertain paraetmers given as 
% 
% $$ k \in [1,10], ~~  h \in [0,1].$$

%%
% This time the numerator has a gain |k| and a delay |h|, and the
% denumerator has a certain double integrator.
% The denumerator thus requires a |poly| qrff member, which encapsulates the 
% certain polynomial $a(s)=s^2$.   
% 

k = qpar('k',5,1,10,8);
h = qpar('h',0,0,1,8);
num = [qrff('gain',k) qrff('delay',h)]; % numerator
den = qrff('poly',[1 0 0]);             % denumerator
P2 = qplant(num,den);
P2.cnom;
P2.ctpl('rff',w);
P2.showtpl

%% 
% One can see that the nominal is drawn at the wrong phase: $+180$ degrees instead
% of $-180$ degrees. This can be corrected using the command |P2.unwrap()|. 
