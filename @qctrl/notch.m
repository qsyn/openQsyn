function obj = notch(wn,wd,zn,zd)

% notch: wn = wd
% skew notch: wn != wd

if nargin==4 
%     zeros_array = roots([1 2*zn*wn wn^2])';
%     poles_array = roots([1 2*zd*wd wd^2])';
%     obj = qctrl(zeros_array, poles_array, 1);
      s = qctrl(0,[],1);
      Czeros = s^2+2*zn*wn*s+wn^2;
      Cpoles = s^2+2*zd*wd*s+wd^2;
      obj = Czeros/Cpoles;
elseif zn > zd
    error('incorrect: should be zn < zd')
else
    error('incorrect number of arguments')    
end