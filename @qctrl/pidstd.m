function obj = pidstd(K,Ti, Td, N)
%pid returns a qctrl pid compensator
%
% Usage:
%
% cpid =  QCTRL.PIDSTD(K, Ti, Td, N) returns a PID compensator in
% standard form
%                 1          Td*s
%   K * ( 1 + ------ + ------------ )
%                Ti*s    (Td/N)*s+1
%
% K: gain
% Ti: Integrator time
% Td: Derivative time
% N: Derivative filter divisor

kp_f = K;
Ti_f = Ti;
Td_f = Td;
N_f = N;

s = qctrl(0,[],1);
%Standard PID
obj = kp_f*(1 + 1/(Ti_f*s) + Td_f*s/((Td_f/N_f)*s+1));
end