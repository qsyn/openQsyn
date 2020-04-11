function [T] = rffpz(obj,w,pzf,dist)

%RFFPZ      produces a real pole or real zero template in  real factored form.
%
%           [T] = rffpz(a,w,form,pzf,dist)
%
%       T: 	column vector with an even number of elements.
%           The upper half the vector T is the low  gain template edge,
%           the second half of the vector T denotes the high  gain
%           template edge (equal in the case of a first order factor).
%           Each element is of the form degree + j*dB
%           (Nichols chart representation). Each edge is sorted
%            in ascending order  w r t angles, in the
%           interval (-180, 180] deg. In some cases the angles are
%           not contiguous (e.g. when the parameter a defines both
%           positive and negative numbers).
%
%       w: 	frequency [rad/s], non-negative real number.
%
%       form: 	The form of the factor, 'dc' or 'hf'.
%               'dc': Indicates dc form (1 + s/a) or 1/(1 + s/a) (default).
%               'hf': Indicates high frequency form (s+a) or 1/(s+a)
%
%
%           pzf: Pole zero flag, 'z' or 'p'.
%               'z': Indicates zero factor (1 + s/a) or (s+a) (default).
%               'p': Indicates pole factor 1/(1 + s/a) or 1/(s+a).
%
%
%       dist: 	The template edges are computed for angles that are
%           multiples of dist [deg] which must be  a positive real
%           number  such that 360 may be divided by dist without
%           remainder The default of dist is 1 [deg]. The angular
%           distance between neigbouring edge points is a multiple of
%           dist.
%
%       End point phase rounding is performed as follows: If a true
%       template end point is nearer a phase grid point outside the
%       the true template, that grid point is included in the computed
%       template, with a gain value equal to that of the true template.
%
%       Reference: Gutman, P-O, Baril C, Neumann L:"An algorithm for
%       computing value sets of uncertain transfer functions in
%       factored real form." IEEE Transactions on Automatic Control,
%       vol 29, no 6, 1268-1273, June 1994.
%
%       See also RFFCPZ, RFFEL, RFFMUL

% Adapted from old Qsyn, Original Authors:
% Author: B Cohen, P-O Gutman
% Version upgrade: A. & Y. Greenhut

if isa(obj.par1,'qpar')
    a = [obj.par1.lower obj.par1.upper]; 
else
    a = [obj.par1 obj.par1];
end
form = obj.type;

% Defaults
if ~(exist('form')==1) ,   form=[];  end
if ~(exist('pzf')==1)  ,   pzf=[];   end
if ~(exist('dist')==1) ,   dist=1;   end

if (abs(360-fix(360/dist)*dist)>eps)
    error('360 must be an integer multiple of the angular resolution, dist [deg]');
end

if isempty(form) , form='dc'; end
if isempty(pzf) ,  pzf='z';   end
if length(a) == 1,     a=[a a];   end

%rffpz handles negative elements in a, too.

a = [min(a) max(a)];

% Find the minimum and the maximum angle of the template  in (0,180) deg
if (a(2) ~= 0)
    p1=rem(atan(w/(a(2)))*180/pi+180,180);
else
    p1=90;
end
if (a(1) ~= 0)
    p2=rem(atan(w/(a(1)))*180/pi+180,180);
else
    p2=90;
end

phim = min(p1,p2);
phiM = max(p1,p2);

% Round to the desired distance.  %Note: phi in (0,180) deg
phi = [ phim round(phim/dist)*dist : dist : round(phiM/dist)*dist phiM]';
n=length(phi);

wc = w./tan(phi*pi/180);  	%Note: tan(pi/2)= 1.6332e+016
sc = size(wc);
s = 1j*w;


if strcmp(form,'dc')
    T=c2n(1j*(w./wc)+1, 180);	%Note: real(T) in (-90,90] but not nec contiguous
    % Note: in Matlab, j/Inf = NaN+i*NaN
else
    T=c2n(s+wc,180); 		%Note: real(T) in (0,180)
end

% phase rounding
if phi(1) > phi(2)
    T(2) = real(T(2))+1j*imag(T(1));
end
if phi(n) < phi(n-1)
    T(n-1) = real(T(n-1))+1j*imag(T(n));
end
T = T(2:n-1);

% %does the phase jump at 90 deg?  %commented out by peo 960620
%   i90 = find(abs(real(T)-90)<1e-10);	%eps does not work
%   if (i90 ~=[]),
%      if ((a(2)==0) & (a(1)<0)),
%        T(i90) = -90 + j*imag(T(i90));
%      end
%      if ((a(2)>0) & (a(1)<0)),
%       T = [T; (-90+j*imag(T(i90)))];	%Note:phases are not well ordered
%     end
%   end

% factor = pole ? zero
if pzf=='p'
    T = -T;
end

% sort phases
[~,ix]=sort(real(T));
T=T(ix);

T=[T ; T];









