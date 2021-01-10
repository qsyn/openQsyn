%% Intro
% This is an example file for the time-domain analysis routines written.
% There are two main routines, qStep and qLsim, which use three main
% subroutines and five auxiliary functions. The auxiliary functions execute
% some basic algorithms that can be useful in future extensions, which is
% why I wrote them separately and as general as possible.
%
% First I'll demonstrate the basic time simulation, and then the specific
% routines themselves.
%% Define qPlant and controllers as in the OpenQsyn example file
k=qpar('k',2,2,5,8);
a=qpar('a',3,1,3,8);
z=qpar('z',0.6,0.3,0.6,8);
wn=qpar('wn',4,4,8,8);
num = [k*wn*wn k*wn*wn*a];
den = [1 2*z*wn wn*wn];
P = qplant(num,den);
den = [1 2*z*wn wn*wn];
P = qplant(num,den);
w_nom = logspace(-2,2,200);
P.cnom(w_nom);


pars = [1 3; 2 5; 8 4; 0.3 0.3];

s = zpk(0,[],1);
set(s,'DisplayFormat','Frequency')
G = 2.5*(1+s/5)*(1+2*0.6*s/4+s^2/16)/s/(1+s)/(1+s/3.2)/(1+s/26);
Pre = 1/(1+2*0.83*s/3.4+s^2/3.4^2);
C = qctrl(G);
F = qctrl(Pre);
%% Define closed loop
L = series(P,C);               % open loop
S = feedback(L,1);             % closed loop from d to y (sensitivity)
T = series(S,series(L,F));     % closed loop from r to y
%% Randomize parameter matrix
rng(0,'twister');
nCases=5;
a_cases = (a.upper-a.lower).*rand(nCases,1) + a.lower;
k_cases = (k.upper-k.lower).*rand(nCases,1) + k.lower;
wn_cases = (k.upper-wn.lower).*rand(nCases,1) + wn.lower;
z_cases = (k.upper-z.lower).*rand(nCases,1) + z.lower;
pars=[a_cases';k_cases';wn_cases';z_cases'];

%% Step response of analog, rational system
y_T=step(T,'Pars',pars);
%% Results w/ Control System toolbox
t=0:0.1:10; %The default time vector in the routine
figure
hold on
for ii=1:nCases
   num_tf = [k_cases(ii)*wn_cases(ii)*wn_cases(ii) k_cases(ii)*wn_cases(ii)*wn_cases(ii)*a_cases(ii)];
   den_tf = [1 2*z_cases(ii)*wn_cases(ii) wn_cases(ii)*wn_cases(ii)]; 
   P_cst=tf(num_tf,den_tf);
   T_cst=minreal(Pre*feedback(P_cst*G,1));
   step(T_cst,t)
end

%% Response to general signal of analog, rational system
t = 0:0.1:30;
u=sin(t);
y=lsim(T,'Pars',pars,'Time',t,'u',u);
%% Results w/ Control System toolbox
t = 0:0.1:30;
figure
hold on
for ii=1:nCases
   num_tf = [k_cases(ii)*wn_cases(ii)*wn_cases(ii) k_cases(ii)*wn_cases(ii)*wn_cases(ii)*a_cases(ii)];
   den_tf = [1 2*z_cases(ii)*wn_cases(ii) wn_cases(ii)*wn_cases(ii)]; 
   P_cst=tf(num_tf,den_tf);
   T_cst=minreal(Pre*feedback(P_cst*G,1));
   lsim(T_cst,u,t)
end
%% How the routines work
% In a nutshell: for each set of parameters we find the
% numerator/denominator data, manipulate the polynomial coefficients to
% obtain the correct transfer function (dictated by the connections
% property of the qSys class), find a state-space realization and simulate.
%
% Break down:
%   1) Parse the parameters: The order of the qpar objects inside a qplant
%       is not identical to the qpar order in the num/den qpoly objects
%       that define it. The utility function Prase_params(obj,pars) outputs
%       two arrays with the parameters in the correct order for the
%       numerator and denominator.
%   2) Find num/den data: CoeffExtract(obj,pars,h) is a recursive function
%       that extracts the correct num/den coefficient vectors from a qSys
%       object. The base case of a qSys object with no "sub" qSys objects
%       in it is solved by ReducedSolver.
%   3) Polynomial manipulations: Inside ReducedSolver the num/den data is
%       extracted, using discrete time convolution and addition (with zero
%       padding when needed) the equivalent transfer function data is
%       obtained according to the connection propety in the qSys object.
%       Note that the routine is flexiable and can be easily modified for
%       different connections/controller classes.
%   4) Find state-space realization: Local_tf2ss takes as input num/den
%       coefficient vectors and outputs [A,B,C,D] matrices. This is a
%       general utility function that might be useful later on, it works
%       only on proper and SISO transfer functions. Obviously both analog
%       and discrete.
%   5) Simulation: Simulation is done via MATLAB's ODE45 solver, qLsim
%       works on general input signals, but this required a bit of
%       "trickery". Since MATLAB does not allow for fixed step solvers
%       outside of Simulink (to my great surprise), a vector of values such
%       as u=sin(t) can't be used directly because the actual solution time
%       vector is unknown a-priori. To solve this, inside the OdeFun I
%       interpolate on u according to the varying time step. It slows the
%       procedure a bit, but works nice all in all.
%% Dealing with delays
% To deal with delays the idea was to discretize the system so delays would
% become simple poles at the origin, and choose a sampling period fast
% enough (say 20 times faster than the delay, although this can be
% changed). Discretization via ZOH equivalent is stable numerically using
% Van-Loan's matrix exponential, Local_c2d implements this algorithm
% directly. However, here lies the main problem: in order for the rest of
% the routine to work, I need to convert this back to rational functions
% (the fact we are now in the Z domain doesn't change anything).
%
% I wrote a function for that, Local_ss2tf, that is supposed to do that.
% Poles are immidiate, I wrote sszero(A,B,C,D) which calculates
% transmission zeros (it is general, works for square MIMO systems as well)
% but for some reason the static gain is still a bit off when compared with
% resuts with Matlab's c2d. This is illustrated below, note that 
% CoeffExtract uses the local routine 'qMinreal' internally, so the 
% resulting transfer function is minimal.
P.delay=[];
h=0.1;
nu=20;
Ts=h/nu;
z_d=tf('z');
Ld_cst=zpk(series(series(c2d(tf(P),Ts),1/z_d^20),c2d(tf(C),Ts)))

adelay(P,h);
L = series(P,C);               % open loop
[LNum,LDen]=CoeffExtract(L,1,h,Ts);
L_qsys=zpk(tf(LNum,LDen,Ts))
%% Effect of the gain drift
% The gain is slightly off, but other then that the results are very close.
% This could pose a problem since discerete systems are gain limited, but
% the routine is modular so this can be fixed easily if we figure out why 
% it is wrong. AbsorbDelay absorbs non fractional delays into the 
% parameters, and all routines have a case for obj.delay ~=0 which 
% discreteizes all the systems and simulates x(k + 1) = A*x(k) + B*u(k)
% instead of \dot{x}(t)=Ax+Bu (with appropriate matrices). This works
% alright for relatively simple systems, although  the gain drift is
% noticable. One advantage we have is that 't's and 'u's sizes are adjusted
% internally according to the discretization, which CST does not do.
t = 0:0.1:30;
u=sin(t);
y=lsim(L,'Pars',pars,'Time',t,'u',u,'Ts',Ts);
close all

h=0.1;
nu=20;
Ts=h/nu;
z_d=tf('z');
s=tf('s');
t2=0:Ts:max(t);
u2=sin(t2);
for ii=1:nCases
   num_tf = [k_cases(ii)*wn_cases(ii)*wn_cases(ii) k_cases(ii)*wn_cases(ii)*wn_cases(ii)*a_cases(ii)];
   den_tf = [1 2*z_cases(ii)*wn_cases(ii) wn_cases(ii)*wn_cases(ii)]; 
   P_cst=tf(num_tf,den_tf);
   Ld_cst=zpk(series(series(c2d(tf(P_cst),Ts),1/z_d^20),c2d(tf(C),Ts)));
   L_cst=series(P_cst*exp(-h*s),G);
   y_dcst(:,ii)=lsim(Ld_cst,u2,t2);
   y_cst(:,ii)=lsim(L_cst,u2,t2);
   %y_cst(:,ii)=step(L_cst,t2);
   %y_dcst(:,ii)=step(Ld_cst,t2);
end
figure
hold on
subplot(1,3,1)
    plot(t2,y_dcst)
    title('CST Discrete system')
subplot(1,3,2)
    plot(t2,y)
    title('OpenQsys Discrete system')
subplot(1,3,3)
    plot(t2,y_cst)
    title('CST Real system')
%% Compound delays
% The effects of the gain mismatch tend to become more pronounces when the
% delays are internal/compound due to feedback interconnections. This is
% not always the case, see for example the response to the set of
% parameters below. The responses aren't identical, but they give the
% "close enough" behaviour. Note that the default sample time for a given
% delay h is Ts=h/20, nicresp isn't defined for qsys objects so there is no
% easy way that I found to make it a function of closed loop bandwidth.
% However, if/when such method is developed this can be adjusted easily.

pars2=[2.629447372786358   2.811583874151238   1.253973632587012   2.826751712278039   2.264718492450819;
   2.292621214998229   2.835494656601145   3.640644557614952   4.872520506302893   4.894665605597829;
   4.157613081677548   4.970592781760615   4.957166948242945   4.485375648722841   4.800280468888801;
   0.966865791547912   2.282278028343492   4.603956968388616   4.023374448929906   4.809614404046644];


S = feedback(L,1);             % closed loop from d to y (sensitivity)
T = series(S,L);     % closed loop from r to y
t = 0:0.1:30;
clear y y_cst y_dcst
u=sin(t);
y=lsim(T,'Pars',pars2,'Time',t,'u',u,'Ts',Ts);
close all

t2=0:Ts:max(t);
u2=sin(t2);

for ii=1:nCases
   %num_tf = [k_cases(ii)*wn_cases(ii)*wn_cases(ii) k_cases(ii)*wn_cases(ii)*wn_cases(ii)*a_cases(ii)];
   %den_tf = [1 2*z_cases(ii)*wn_cases(ii) wn_cases(ii)*wn_cases(ii)]; 
   num_tf = [pars2(2,ii)*pars2(3,ii)*pars2(3,ii) pars2(2,ii)*pars2(3,ii)*pars2(3,ii)*pars2(1,ii)];
   den_tf = [1 2*pars2(4,ii)*pars2(3,ii) pars2(3,ii)*pars2(3,ii)]; 
   P_cst=tf(num_tf,den_tf);
   Ld_cst=zpk(series(series(c2d(tf(P_cst),Ts),1/z_d^20),c2d(tf(C),Ts)));
   Td_cst=zpk(feedback(Ld_cst,1));
   L_cst=series(P_cst*exp(-h*s),G);
   T_cst=feedback(L_cst,1);
   y_dcst(:,ii)=lsim(Td_cst,u2,t2);
   y_cst(:,ii)=lsim(T_cst,u2,t2);
end
figure
hold on
subplot(1,3,1)
    plot(t2,y_dcst)
    title('CST Discrete system')
subplot(1,3,2)
    plot(t2,y)
    title('OpenQsys Discrete system')
subplot(1,3,3)
    plot(t2,y_cst)
    title('CST Real system')
%% Sampled systems
% qPlant objects don't seem to have a sample-time property, and I didn't
% address simulating hybrid systems directly (via discretizing everything),
% but the current routines in place can be adapted to it. They already
% recieve sample time as an external parameter, and Local_c2d does ZOH
% equivalent transformations.