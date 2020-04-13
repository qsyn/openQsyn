function [T] = rffel(obj,w,dist)

%RFFEL      pure gain, delay, unstructured uncertainty, or integrators rff template
%
%   [T] = rffel(element,a,w,dist)
%
%   produces a
%   template in  real factored form for uncertain gain,
%   uncertain delay, multiplicative unstructured uncertainty,
%   and an uncertain number of integrators.
%   T: column vector with an even number of elements.
%   The first (upper, with lowest indeces)  half the vector T denotes
%   the low gain template border,
%   the second half of the vector T denotes the low gain
%   template border. Each vector element is of the form degree + j*dB
%   (Nichols chart representation).
%
%   element: 'gain', 'delay', 'uns', or 'int'. See below.
%
%   a: Vector with two or more real elements. See below.
%
%   w: frequency [rad/s], non-negative real number, for which the
%   template is computed.
%
%   dist: The phases of the computed template values are integer
%   multiples of dist [deg]. The angular distance between two
%   neighbouring computed template point is dist [deg] (default = 1).
%   dist must be such that 360 is divisible by dist without remainder
%
%[T] = rffel('gain',a,w,dist) produces a magnitude template.
%
%   a: Vector with two real elements in the form of
%   a=[amin amax]. amin represents the minimum gain [absolute
%   value] and amax represents the maximum gain [absolute value].
%   Maximum and minimum gain must have the same non-zero sign.
%
%[T] = rffel('delay',a,w,dist) produces a delay template.
%
%   a: Vector with two non-negative real numbers in the form of
%   a=[amin amax]. amin represents the minimum delay [seconds]
%   and amax represents the maximum delay [seconds].
%
%[T] = rffel('uns',m,dist,w) produces an unstructured multiplicative
%   uncertainty template.
%
%   m: Matrix with two rows,   or  alternatively,  one real number in [0,1) .
%       If m is a matrix, then m must contain at least two columns.
%       The first row, containing non-negative real numbers holds
%       frequencies [rad/s] in increasing order. The second row,
%       consisting of real numbers in the interval [0, 1),
%       contains the unstructured multiplicative uncertainty radii
%       [absolute value], for each of the frequencies in the first
%       row, respectively. The uncertainty radius is then linearly
%       interpolated, or end point constant extrapolated, respectively,
%       with respect to frequency.
%       If m is a single number, then it denotes the uncertainty
%       radius for all frequencies
%
%[T] = rffel('int',a,w,dist) produces a templete for an uncertain number
%   of differentiators/integrators.
%
%   a: Vector of two integer numbers, in the form of
%   a=[amin amax]. amin represents the minimum number of integrators,
%   and amax represents the maximum  maximum number of integrators.
%   (-amin represents the maximum number of differentiators, and
%   -amax represents the minimum number of differentiators)
%
%   End point phase rounding is performed as follows: If a true
%   template end point is nearer a phase grid point outside the
%   the true template, that grid point is included in the computed
%   template, with a gain value equal to that of the true template.
%
%   Reference: Gutman, P-O, Baril C, Neumann L:"An algorithm for
%   computing value sets of uncertain transfer functions in
%   factored real form." IEEE Transactions on Automatic Control,
%   vol 29, no 6, 1268-1273, June 1994.
%
%
%
%
%   See also RFFPZ, RFFCPZ, RFFMUL

% Adapted from old Qsyn, Original Authors:
% Authors:	B. Cohen, P-O Gutman
% Version upgrade: A. & Y. Greenhut

if isa(obj.par1,'qpar')
    a = [obj.par1.lower obj.par1.upper]; 
else
    a = [obj.par1 obj.par1];
end
element = obj.type;

if ~exist('dist','var'), dist=1; end

if strcmp(element,'gain')              % Gain
    
    t=[min(a) max(a) ]';
    if sign(t(1))~=sign(t(2)) || t(1)==0	% sign(0)=0
        error(' max and min gain must have the same non-zero sign')
    end
    
    if sign(t(1))>0
        phi = 0;
    else
        phi = round(-180/dist)*dist;
    end
    T=[phi;phi] + 20*log10(abs(t))*1j;
end

if strcmp(element,'delay')           % t=e^(-a*s)
 
    phim = -(max(a)*w*180/pi);		%least phase, degrees
    phiM = -(min(a)*w*180/pi);
    T = [round(phim/dist)*dist : dist : round(phiM/dist)*dist]'; %column vector
    T=[T ; T];
end

if strcmp(element,'uns')            % t=1+M(w)

    if length(a) == 1
        if a<0 || a>=1
            error(' modulus of unstructured uncertainty must be in [0,1)');
        end
        M = a;
    else
        if min(a(1,:))<0
            error(' frequencies for unstructured uncertainty must be positive!');
        end
        if min(a(2,:))<0 || max(a(2,:))>=1
            error(' modulus of unstructured uncertainty must be in [0,1)');
        end
        if  min(size(a))<2
            if a(1,1)~=w
                error('unstructured uncertainty description matrix must have at least two columns')
            else
                a = [a;a];
            end
        end
        if sort(a(1,:))~=a(1,:)
            error(' frequencies for unstructured uncertainty is not given in ascending order');
        end
        if a(1,1)>w
            M = a(2,1);
        elseif a(1,length(a(1,:)))<w
            M = a(2,length(a(1,:)));
        else
            M=interp1(a(1,:),a(2,:),w,'linear');
        end
    end
    phim = -asin(M)*180/pi;
    phiM =  asin(M)*180/pi;
    
    phi = [round(phim/dist)*dist : dist : round(phiM/dist)*dist];
    phi = phi(find( (phi-phim>=0)&(phi-phiM<=0) )); % desired angles inside the interval
    
    phi = [ phim phi phiM]';	%contains at least three elements, since 0 is in
    
    
    v = sqrt(M*M-sin(phi*pi/180).*sin(phi*pi/180));
    
    r1=cos(phi*pi/180)-v;	%distances from origin to intersection with unc circle
    r2=cos(phi*pi/180)+v;
    
    %phase rounding
    n=length(phi);
    if phi(1) < (phi(2) - dist/2) %note: symmetric angles
        phi(n) = phi(n-1) + dist;
        phi(1) = phi(2)-dist;
    else
        phi(n)=[]; r1(n)=[]; r2(n)=[];
        phi(1)=[]; r1(1)=[]; r2(1)=[];
    end
    
    
    r=[r1 ; r2];
    %t=r.*exp([phi; phi ]*pi/180*j);
    T= [phi;phi] + 20*log10(abs(r))*1j;
    
end

if strcmp(element,'int')          % t=1/(s^n)  n== integer;
    
    if max(abs(rem(a,1)))>1e-10
        error('the number of integrators or differentiators is not integer')
    end
    s=1j*w;
    n=[min(a):1:max(a)]' ;
    
    T= round(-90*n/dist)*dist + 1j*20*log10(abs((1)./s.^n));
    T=[T;T];
    
end



