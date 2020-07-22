function obj = pid(kp, ki, kd, Tf)
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
switch nargin
    case 1
        kp_f = kp;
        ki_f = 0;
        kd_f = 0;
        Tf_f = 0;
    case  2
        kp_f = kp;
        ki_f = ki;
        kd_f = 0;
        Tf_f = 0;
    case 3
        kp_f = kp;
        ki_f = ki;
        kd_f = kd;
        Tf_f = 0;
    case 4
        kp_f = kp;
        ki_f = ki;
        kd_f = kd;
        Tf_f = Tf;    
end
s = qctrl(0,[],1);
%PID
    obj = kp_f + ki_f/s + (kd_f*s)/(Tf_f*s +1); % with 1st order derivative filter
end
