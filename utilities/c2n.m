function [y]=c2n(x,op)

%C2N        converts a complex matrix from complex to Nichols form.
%           y = c2n(x,op) 
%
%   Outputs:
%
%   y:  resulting matrix in Nichols form, i.e. each element is given as
%       deg + j*dB. 
%   x:  original matrix in complex form, i.e. each element is given as
%       a + j*b. 
%   op: one of the following:
%           'wrap' = for each matrix element, the phase, deg, is 
%               wrapped into the range [-360,0] degrees (default);
%           'unwrap' = the phase is unwrapped columnwise, i.e. 
%               continuous over Riemann surface boundaries, 
%               if possible;	
%           real number, c = the phase is unwrapped columnwise, i.e. 
%               continuous over Riemann surface boundaires, if 
%               possible, with the phase of the (1,1)-matrix element 
%               belonging to the Riemann surface [c-180,c+180] degrees.
%
%   Example:
%   a = [
%           1.0000 + 0.1000i   1.0000 - 0.1000i
%           1.0000 - 0.1000i  -1.0000 - 0.1000i
%           -1.0000 - 0.1000i  -1.0000 + 0.1000i
%           -1.0000 + 0.1000i   1.0000 + 0.1000i];
%   y=c2n(a)
%
%   y =
%           1.0e+002 *
%           -3.5429 + 0.0004i  -0.0571 + 0.0004i
%           -0.0571 + 0.0004i  -1.7429 + 0.0004i
%           -1.7429 + 0.0004i  -1.8571 + 0.0004i
%           -1.8571 + 0.0004i  -3.5429 + 0.0004i
%   y=c2n(a,'unwrap')
%
%   y =
%           1.0e+002 *
%           0.0571 + 0.0004i  -0.0571 + 0.0004i
%           -0.0571 + 0.0004i  -1.7429 + 0.0004i
%           -1.7429 + 0.0004i  -1.8571 + 0.0004i
%           -1.8571 + 0.0004i  -3.5429 + 0.0004i
%
%  
%   y=c2n(a,-179)
%
%   y =
%           1.0e+002 *
%           -3.5429 + 0.0004i  -3.6571 + 0.0004i
%           -3.6571 + 0.0004i  -5.3429 + 0.0004i
%           -5.3429 + 0.0004i  -5.4571 + 0.0004i
%           -5.4571 + 0.0004i  -7.1429 + 0.0004i
%
%
%
%   See also: Matlab function UNWRAP, N2C
%
%
%


% Author: B Cohen, M Nordin, P-O Gutman 
% Copyright: El-Op Electro-Optics Industries Ltd, and Novomatic AB 1996
% Version upgrade: A. & Y. Greenhut


if nargin==0 ,
   disp('  [y]=c2n(x,op)')
   return;
end;

y=[];
if isempty(x);return;end;																						

% Defaults
if ~(exist('op')==1), op='wrap'; end; %J

if strcmp(op,'wrap'),
   y=rem(angle(x)*180*(1/pi)+10*360,360)-360+20*log10(abs(x))*j;
   return;
else,
   index1 = find(~isnan(x));
   %ax1 = unwrap(angle(x(index1)));
   ax1 = unwrap(angle(x(index1)),[],2);
   ax = zeros(size(x));
   ax(index1) = ax1;

   y=ax*(180/pi)+20*log10(abs(x))*1j;
   if ~isstr(op),
      n_r=ceil((real(y(1,1))-(op+180))/360);
      y=y-(n_r)*360;
   end;
end;

        
        


