function obj = pid(kp, ki, kd, Tf, Ti, Td, N, flag)
%pid returns a qctrl pid compensator
%
% Usage:
%
% cpid = qctrl.pid(kp, ki, kd, Tf, Ti, Td, N)  returns a pid compensator,
% based on the input values. there are two cases: pid with kp,ki,kd and Tf
% or standard pid with Ti,Td and N. It is allow to create only PI for
% example, with qctrl.pid(kp,ki).
% the inputs are:
% kp: Proportional gain
% ki: Integral gain
% kd: Derivative gain
% Tf: Time constant of the first-order derivative filter
% Ti: Integrator time
% Td: Derivative time
% N: Derivative filter divisor
% flag: selector for pid or standard pid

switch nargin
    case 1
        f_flag = 0;
        kp_f = kp;
        ki_f = 0;
        kd_f = 0;
        Tf_f = 0;
        Ti_f = 0;
        Td_f = 0;
        N_f = 0;
    case  2
        f_flag = 0;
        kp_f = kp;
        ki_f = ki;
        kd_f = 0;
        Tf_f = 0;
        Ti_f = 0;
        Td_f = 0;
        N_f = 0;
    case 3
        f_flag = 0;
        kp_f = kp;
        ki_f = ki;
        kd_f = kd;
        Tf_f = 0;
        Ti_f = 0;
        Td_f = 0;
        N_f = 0;
    case 4
        f_flag = 0;
        kp_f = kp;
        ki_f = ki;
        kd_f = kd;
        Tf_f = Tf;
        Ti_f = 0;
        Td_f = 0;
        N_f = 0;
    case 8
        f_flag = flag;
        kp_f = kp;
        ki_f = ki;
        kd_f = kd;
        Tf_f = Tf;
        Ti_f = Ti;
        Td_f = Td;
        N_f = N;        
end

s = qctrl(0,[],1);

%PID

if f_flag == 0
    obj = kp_f + ki_f/s + (kd_f*s)/(Tf_f*s +1); % with 1st order derivative filter
end

%Standard PID

if f_flag == 1
    obj = kp_f*(1 + 1/(Ti_f*s) + Td_f*s/((Td_f/N_f)*s+1));
end
if f_flag ~=1 && f_flag ~=0
    error('not enough arguments')
end
end