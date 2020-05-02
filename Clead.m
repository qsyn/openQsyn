function Clead = Clead(DesiredAngle,Freq,Order)
%this function output the Clead, when its inputs are the desired angle,
%frequency and order
s = tf('s');
wn = Freq;
syms alpha1
eqn = DesiredAngle == asind((alpha1 - 1)/(alpha1 + 1));
a = double(solve(eqn,alpha1));
switch Order
    case {1}
        Clead = (sqrt(a)*s+wn)/(s+sqrt(a)*wn);
    case {2}
       Clead = (sqrt(a)*s+wn)/(s+sqrt(a)*wn); % still need to fix this 
end
end