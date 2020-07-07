function obj = notch(wn,wd,zn,zd)
%notch returns a QCTRL notch compensator
%
% Usage:
%
% cnotch = qctrl.notch(wn,wd,zn,zd)  returns a notch compensator
% at a given frequency.
%
% notch:zd > zn
% anti-notch : zn>zd
% skew notch: wn != wd

if nargin==4
    zero_1  = -zn*wn + sqrt((wn^2)*(zn-1));
    zero_2  = -zn*wn - sqrt((wn^2)*(zn-1));
    pole_1  = -zd*wd + sqrt((wd^2)*(zd-1));
    pole_2  = -zd*wd - sqrt((wd^2)*(zd-1));
    zeros = [zero_1 zero_2];
    poles = [pole_1 pole_2];
        obj = qctrl(zeros,poles,1);
else
    error('incorrect number of arguments')
end