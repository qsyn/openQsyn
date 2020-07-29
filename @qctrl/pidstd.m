<<<<<<< HEAD
function obj = pidstd(K,Ti, Td, N)
=======
function obj = pidstd(kp,Ti, Td, N)
>>>>>>> tests
%pid returns a qctrl pid compensator
%
% Usage:
%
<<<<<<< HEAD
% cpid =  QCTRL.PIDSTD(K, Ti, Td, N) returns a PID compensator in
% standard form
%                 1          Td*s
%   K * ( 1 + ------ + ------------ )
%                Ti*s    (Td/N)*s+1
%
% K: gain
=======
% cpid = qctrl.pidstd(Ti, Td, N)  returns a standard pid compensator,
% based on the input values.
% the inputs are:
>>>>>>> tests
% Ti: Integrator time
% Td: Derivative time
% N: Derivative filter divisor

<<<<<<< HEAD
kp_f = K;
=======
kp_f = kp;
>>>>>>> tests
Ti_f = Ti;
Td_f = Td;
N_f = N;

s = qctrl(0,[],1);
%Standard PID
obj = kp_f*(1 + 1/(Ti_f*s) + Td_f*s/((Td_f/N_f)*s+1));
end