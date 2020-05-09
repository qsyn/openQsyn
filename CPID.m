function CPID = CPID(ki, kp, kd, Tf, taud, flag) 

% kp: Proportional gain
% ki: Integral gain
% kd: Derivative gain
% Tf: Time constant of the first-order derivative filter
% taui: Reset time
% taud: Derivative time
% flag: 0 - gain design , 1 - time design, 2 - Ziegler-Nichols tuning

s = tf('s');
if flag == 0
    if nargin <= 4 
        Tf = 0; % the controller has no filter on the derivative action 
    end
    CPID = kp + ki/s + (kd*s)/(Tf*s +1);
end

if flag == 1
    taui = 1/ki;
    CPID = kp*(1 + 1/(taui*s) + taud*s);
end

% if flag == 2
%     % find a way to change kp with respect to the amplite of the system's
%     % output, in such a way that average of the amplited is zero. Then,
%     % apply the Z-N tuning rules with ku,tauu
% end

end