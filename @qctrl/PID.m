function obj = PID(kp, ki, kd, Tf, Ti, Td, N, flag) 

% kp: Proportional gain
% ki: Integral gain
% kd: Derivative gain
% Tf: Time constant of the first-order derivative filter
% Ti: Integrator time
% Td: Derivative time
% N: Derivative filter divisor
% flag: 0 - PID , 1 - Standard PID

s = qctrl(0,[],1);

%PID

if flag == 0
    obj = kp + ki/s + (kd*s)/(Tf*s +1); % with 1st order derivative filter
end

%Standard PID

if flag == 1
    obj = kp*(1 + 1/(Ti*s) + Td*s/((Td/N)*s+1));
end
if flag ~=1 && flag ~=0
    error('not enough arguments')
end
end