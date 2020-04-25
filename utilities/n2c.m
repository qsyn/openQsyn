function[y]=n2c(z);

%N2C        converts a matrix from Nichols form to complex form
%
%               [y]=n2c(z);
%	
%       y:      resulting matrix in complex form, i.e. each element is given as
%           a + j*b. 
%
%       x:      original matrix in Nichols form, i.e. each element is given as
%               deg + j*dB. 


% This file is part of the OpenQsyn toolox, distribted under GNU LGPL3 license
% Author: B. Cohen
% Version Upgrade: A. & Y. Greenhut

y=10.^(0.05*imag(z)).*exp(j*real(z)*(pi/180));

