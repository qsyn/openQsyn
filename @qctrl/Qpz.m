function obj = Qpz(a,b,c,d,flag)
s = qctrl(0,[],1);
obj = 1;
switch flag
    case 0 
        obj = 1/s;
    case 1
        obj = s;
    case 2 
        obj = 1/(s+a);
    case 3
        obj = 1/(s^2+a);
    case 4
        obj = s^2+a;
    case 5
        obj = s+a;
    case 6
        obj = (c*s+d)/(a*s+b);
end
