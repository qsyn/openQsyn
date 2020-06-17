function Cpz = Cpz(a,b,c,d,flag)
s = qctrl(0,[],1);
Cpz = 1;
switch flag
    case 0 
        Cpz = 1/s;
    case 1
        Cpz = s;
    case 2 
        Cpz = 1/(s+a);
    case 3
        Cpz = 1/(s^2+a);
    case 4
        Cpz = s^2+a;
    case 5
        Cpz = s+a;
    case 6
        Cpz = (c*s+d)/(a*s+b);
end
