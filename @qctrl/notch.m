function obj = notch(wn,wd,zn,zd)

% notch: wn = wd
% skew notch: wn != wd

if nargin==4
    %     zeros_array = roots([1 2*zn*wn wn^2])';
    %     poles_array = roots([1 2*zd*wd wd^2])';
    %     obj = qctrl(zeros_array, poles_array, 1);
    %       s = qctrl(0,[],1);
    %       Czeros = s^2+2*zn*wn*s+wn^2;
    %       Cpoles = s^2+2*zd*wd*s+wd^2;
    zero_1  = -zn*wn + sqrt((wn^2)*(zn-1));
    zero_2  = -zn*wn - sqrt((wn^2)*(zn-1));
    pole_1  = -zd*wd + sqrt((wd^2)*(zd-1));
    pole_2  = -zd*wd - sqrt((wd^2)*(zd-1));
    zeros = [zero_1 zero_2];
    poles = [pole_1 pole_2];
%     obj = Czeros/Cpoles;
        obj = qctrl(zeros,poles,1);
else
    error('incorrect number of arguments')
end