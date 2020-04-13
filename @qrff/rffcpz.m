function [T]=rffcpz(obj,w,pzf,angDist)

%RFFCPZ     computes a complex pole/zero pair template in real factored form.
%
%       [T] = RFFCPZ(a1,a2,w,form,pzf,angDist)
%
%       T: 	column vector with an even number of elements.
%       	The first half the vector T is the low gain template edge,
%           the second half of the vector T is the high gain template edge.
%       	Each element is of the form degree + j*dB. Each template edge
%       	is sorted in ascending order with respect to angle.
%       	The same angles occur in both edges.  T contains
%       	angles in the interval (-180, 180] deg, but may have angles
%       	outside this interval in order to make the angle sequence
%       	contiguous. Each angle is a multiple of dist [deg], see below.
%
%       a1: 	When  uncertain, a1 is a vector with two real elements
%       	in the form a1=[zmin zmax]. zmin and zmax represent the
%       	minimum and maximum values of the uncertain relative
%       	damping z. When the relative damping is certain,
%       	either zmin = zmax, or a1 contains one element only.
%       	The elements of a1 are real numbers.
%
%       a2: 	When  uncertain, a2 is a vector with two positive real
%       	elements in the form a2=[wnmin wnmax]. wnmin and wnmax
%       	represent the minimum and maximum values of the uncertain
%       	resonance frequency wn [rad/s]. When wn is certain, either
%       	wmin = wmax, or a2 contains one element only.
%       	The elements of a2 are positive real numbers.
%
%       w: 	[rad/s], non-negative real number, the freqeuncy for which the
%           is computed.
%
%       form: 	The form of the factor, 'dc' or 'hf.
%               'dc': 	Indicates dc form (1 + 2*z*s/wn + s^2/wn^2)
%               or 1/(1 + 2*z*s/wn + s^2/wn^2) (default).
%               'hf': 	Indicates high frequency form (s^2 + 2*z*wn*s + wn^2)
%               or 1/(s^2 + 2*z*wn*s + wn^2).
%
%		pzf:    Pole/zero flag, 'z' or 'p'.
%               'z': 	Indicates zero factor (1 + 2*z*s/wn + s^2/wn^2)
%               or (s^2 + 2*z*wn*s + wn^2) (default).
%               'p': 	Indicates pole factor 1/(1 + 2*z*s/wn + s^2/wn^2)
%               or 1/(s^2 + 2*z*wn*s + wn^2).
%
%
%       angDist:   The template edges are computed for angles that are
%       	multiples of dist [deg] which must be  a positive real
%       	number  such that 360 may be divided by dist without
%       	remainder The default of dist is 1 [deg]. The angular
%       	distance between neigbouring edge points is a multiple of
%       	dist.
%
%       End point phase rounding is performed as follows: If a true
%       template end point is nearer a phase grid point outside the
%       the true template, that grid point is included in the computed
%       template, with a gain value equal to that of the true template end point.
%
%       Reference: Gutman, P-O, Baril C, Neumann L: "An algorithm for
%       computing value sets of uncertain transfer functions in
%       factored real form." IEEE Transactions on Automatic Control,
%       vol 29, no 6, 1268-1273, June 1994.
%
%       See also RFFPZ, RFFEL, RFFMUL

% Adapted from old Qsyn,
% Original Authors: P-O Gutman, B Cohen
% Version update: A. & Y. Greenhut

if ~exist('angDist','var'), angDist=1; end
if isa(obj.par1,'qpar')
    a2 = [obj.par1.lower obj.par1.upper]; % wn
else
    a2 = [obj.par1 obj.par1];
end
if isa(obj.par2,'qpar')
    a1 = [obj.par2.lower obj.par2.upper]; % zeta
else
    a1 = [obj.par2 obj.par2];
end
form = obj.type;

if (abs(360-fix(360/angDist)*angDist)>eps)
    error('360 must be an integer multiple of the angular resolution, dist [deg]');
end

if ( (isempty(a1)) || (isempty(a2)) || (length(a1)>2) || (length(a2)>2) )
    error('The number of relative damping and natural frequency parameters is wrong');
end

if (min(a2)<=0)
    error('Natural frequencies are positive!');
end

% Defaults
if ~exist('form','var'), form=[]; end
if ~exist('pzf','var'),  pzf=[];  end
if ~exist('angDist','var'), angDist=1;  end

if isempty(form),form='dc'; end
if isempty(pzf), pzf='z';   end
if length(a1)==2, a1 = [min(a1) max(a1)]; end
if length(a2)==2, a2 = [min(a2) max(a2)]; end
if length(a1)==1
    a=[a1 a1];
else
    a = a1;
end
if length(a2)==1
    a=[a a2 a2];
else
    a=[a a2];
end


s=1j*w;

% Check if it is a certain case
% =============================

if length(a1)==1, a1(2)=a1(1); end	% peo 960518
if length(a2)==1, a2(2)=a2(1); end	% peo 960518
if (  ( (a1(1)==a1(2)) && (a2(1)==a2(2)))) % peo 960518
    % if ( ( (length(a1)==1) & (length(a2)==1) ) | ( (a1(1)==a1(2)) & (a2(1)==a2(2)))),
    
    wn=a(3); z=a(1);
    if strcmp(form,'dc')
        t=((s*s)/(wn*wn)+2*(z/wn)*s+1);
    else
        t=(s*s+2*z*wn*s+wn*wn);
    end
    
    if ((abs(wn-w)<=eps) && (abs(z)<=eps)) 	%avoid matlab's random phase
        T = c2n(eps);
    else
        T=c2n(t);
    end
    
    if strcmp(pzf,'p')
        T=-T;
    end
    if (real(T)>180), T=T-360; end
    if (real(T)<=-180), T=T+360; end   %Since T is certain, T need not have
    % an upper an lower edge. see rffmul
    
    %round it again, man, because c2n may introduce a slight error:
    T = round(real(T)/angDist)*angDist + 1j*imag(T);   %peo 961015
    
    return
end

% a = [zmin zmax wmin wmax] = [a1 a2]
zmin=a(1);
zmax=a(2);
wmin=a(3);
wmax=a(4);
z=[];
wn=[];

% case 1: (zmin,wmin)-->(zmin,wmax)
% ________________________________
phi_=angle(s*s+2*zmin*[wmin wmax]*s+[wmin wmax].^2)*180/pi; %note: angle(-1)=pi
if ((zmin==0) && (wmin==w)), phi_ = 0; end
if ((zmin==0) && (wmax==w)), phi_ = 180; end
phi11=min(phi_);
phi12=max(phi_);
phi = [round(phi11/angDist)*angDist : angDist : round(phi12/angDist)*angDist];
phi = phi(find( (phi-phi11>=0)&(phi-phi12<=0) )); % desired angles inside the interval

% calculating extreme points
Tmax1 = qrff.rffutil1(w,[phi11 phi12],zmin,zmax,wmin,wmax,form,pzf,1);
% calculate template points
if (~isempty(phi))
    T1 = qrff.rffutil1(w,phi,zmin,zmax,wmin,wmax,form,pzf,1);
else
    T1=[];
end


% case 2: (zmin,wmax)-->(zmax,wmax)
% _________________________________
phi_=angle(s^2+2*[zmin zmax]*wmax*s+wmax^2)*180/pi;
if (wmax==w)
    phi_=angle(2*[zmin zmax]*wmax*s)*180/pi;
    if (zmin==0), phi_(1)=[]; end
    if (zmax==0), phi_(2)=[]; end
end
phi21=min(phi_);
phi22=max(phi_);
if ( ((wmax-w)<0) && (zmin<0) && (zmax>=0) )	%angle(-1)=180 deg
    phi21 = max(phi_);			%minimum angle
    phi22 = min(phi_) + 360;		%maximum angle
end
phi = [round(phi21/angDist)*angDist : angDist : round(phi22/angDist)*angDist];
phi = phi(find( (phi-phi21>=0)&(phi-phi22<=0) )); % desired angles inside the interval

% calculating extreme points
Tmax2 = qrff.rffutil1(w,[phi21 phi22],zmin,zmax,wmin,wmax,form,pzf,2);
% calculate template points
if (~isempty(phi))
    T2 = qrff.rffutil1(w,phi,zmin,zmax,wmin,wmax,form,pzf,2);
else
    T2=[];
end


% case 3: (zmax,wmax)-->(zmax,wmin)
% _________________________________
phi_=angle(s*s+2*zmax*[wmin wmax]*s+[wmin wmax].^2)*180/pi;
if ( (zmax==0) && (wmin==w) ), phi_ = 0; end
if ( (zmax==0) && (wmax==w) ), phi_ = 180; end
phi31=min(phi_);
phi32=max(phi_);
phi = [round(phi31/angDist)*angDist : angDist : round(phi32/angDist)*angDist];
phi = phi(find( (phi-phi31>=0)&(phi-phi32<=0) )); % desired angles inside the interval

% calculating extreme points
Tmax3 = qrff.rffutil1(w,[phi31 phi32],zmin,zmax,wmin,wmax,form,pzf,3);
% calculate template points
if (~isempty(phi))
    T3 = qrff.rffutil1(w,phi,zmin,zmax,wmin,wmax,form,pzf,3);
else
    T3=[];
end



% case 4: (zmax,wmin)-->(zmin,wmin)
% _________________________________
phi_=angle(s^2+2*[zmax zmin]*wmin*s+wmin^2)*180/pi;
if (wmin==w)
    phi_=angle(2*[zmin zmax]*wmin*s)*180/pi;
    if (zmin==0), phi_(1)=[]; end
    if (zmax==0), phi_(2)=[]; end
end
phi41=min(phi_);
phi42=max(phi_);
if ( ((wmin-w)<0) && (zmin<0) && (zmax>=0) )	%angle(-1)=180 deg
    phi41 = max(phi_);			%minimum angle
    phi42 = min(phi_) + 360;		%maximum angle
end
phi = [round(phi41/angDist)*angDist : angDist : round(phi42/angDist)*angDist];
phi = phi(find( (phi-phi41>=0)&(phi-phi42<=0) )); % desired angles inside the interval

% calculating extreme points
Tmax4 = qrff.rffutil1(w,[phi41 phi42],zmin,zmax,wmin,wmax,form,pzf,4);
% calculate template points
if (~isempty(phi))
    T4 = qrff.rffutil1(w,phi,zmin,zmax,wmin,wmax,form,pzf,4);
else
    T4=[];
end



% _________________________________
% Note that rffutil1 gave   template segments for 'z' and in (-180,180] deg
% ---------------------------------


% The origin belongs to the interior of the template
% _________________________________
if ( (sign((w-wmin)*(wmax-w)) > 0) && (sign(zmin*zmax) < 0) )
    Thi = [ T1(:) ; T2(:) ; T3(:) ; T4(:) ];	%column vector
    %disp(' test 1234kv0')
    Thi = rffutil3([],Thi,[],'hi',angDist);	%sort and eliminate angular doubles
    
    if strcmp(form,'hf')		%get gain right for pole-zero cancellations
        Tlo = real(Thi) + 1j*imag(c2n(eps));
    else
        Tlo = real(Thi) + 1j*imag(c2n(eps/(w*w)));
    end
    
    
    % The template lies in the first and second quadrants, with the origin on one edge
    % _________________________________
elseif ( (sign((w-wmin)*(wmax-w)) > 0) && ( (zmin>= 0) && (zmin < eps) ) )
    Thi = [T2(:) ; T3(:) ; T4(:) ];	%high gain edge column vector
    %disp(' test 12kv0')
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = rffutil3(Tmax2,Thi,Tmax4,'hi',angDist);
    
    if strcmp(form,'hf')		%get gain right for pole-zero cancellations
        Tlo = real(Thi) + 1j*imag(c2n(eps));
    else
        Tlo = real(Thi) + 1j*imag(c2n(eps/(w*w)));
    end
    
    
    % The template lies in the third and fourth quadrants, with the origin on one edge
    % _________________________________
elseif ( (sign((w-wmin)*(wmax-w)) > 0) && ( (zmax<= 0) && (zmax > -eps) ) )
    Thi = [T4(:) ; T1(:) ; T2(:) ];	%high gain edge column vector
    %disp(' test 34kv0')
    %change +180 deg to -180 deg to get contiguous phases
    ix = find(abs(real(Thi)-180) < 1e-10);
    if (~isempty(ix)), Thi(ix)= -180 + 1j*imag(Thi(ix)); end
    Tmax4(1) = -180 + 1j*imag(Tmax4(1));
    
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = rffutil3(Tmax4,Thi,Tmax2,'hi',angDist);
    
    if strcmp(form,'hf')		%get gain right for pole-zero cancellations
        Tlo = real(Thi) + 1j*imag(c2n(eps));
    else
        Tlo = real(Thi) + 1j*imag(c2n(eps/(w*w)));
    end
    
    
    %plot(c2n(-w*w + 2*j*w*[zmin:.1:zmax]'*[wmin:0.1:wmax] + ones(size([zmin:.1:zmax]'))*([wmin:0.1:wmax].*[wmin:0.1:wmax]),0),'og')
    
    % The template lies in the fourth and first quadrants, with the origin on one edge
    % _________________________________
elseif ( ( (wmin-w>=0) && ((wmin-w)<eps) ) && (sign(zmin*zmax)<0) )
    Thi = [T1(:) ; T2(:) ; T3(:) ];	%high gain edge column vector
    %disp(' test 41kv0')
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = rffutil3(Tmax1,Thi,Tmax3,'hi',angDist);
    
    if strcmp(form,'hf')	%get gain right for pole-zero cancellations
        Tlo = real(Thi) + 1j*imag(c2n(eps));
    else
        Tlo = real(Thi) + 1j*imag(c2n(eps/(w*w)));
    end
    
    
    % The template lies in the second and third quadrants, with the origin on one edge
    % _________________________________
elseif ( ( (wmax-w<=0) && ((w-wmax)<eps) ) && (sign(zmin*zmax)<0) )
    Thi = [T3(:) ; T4(:) ; T1(:) ];	%high gain edge column vector
    %disp(' test 23kv0')
    %add +360 deg to negative phases to get contiguous phases
    ix = find(real(Thi) < 0);
    if (~isempty(ix)), Thi(ix)= Thi(ix) + 360; end
    Tmax1  = Tmax1 + 360;
    
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = rffutil3(Tmax3,Thi,Tmax1,'hi',angDist);
    
    if strcmp(form,'hf')		%get gain right for pole-zero cancellations
        Tlo = real(Thi) + 1j*imag(c2n(eps));
    else
        Tlo = real(Thi) + 1j*imag(c2n(eps/(w*w)));
    end
    
    
    % The template lies in the  first quadrant, with the origin on one vertex
    % _________________________________
elseif ( ( (wmin-w>=0) && ((wmin-w)<eps) ) && ( (zmin>= 0) && (zmin < eps) ) )
    Thi = [ T2(:) ; T3(:) ];	%high gain edge column vector
    %disp(' test 1kv0')
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = rffutil3(Tmax2,Thi,Tmax3,'hi',angDist);	%
    
    if strcmp(form,'hf')		%get gain right for pole-zero cancellations
        Tlo = real(Thi) + 1j*imag(c2n(eps));
    else
        Tlo = real(Thi) + 1j*imag(c2n(eps/(w*w)));
    end
    
    
    % The template lies in the  second quadrant, with the origin on one vertex
    % _________________________________
elseif ( ( (w-wmax>=0) && ((w-wmax)<eps) ) && ( (zmin>= 0) && (zmin < eps) ) )
    Thi = [ T3(:) ; T4(:) ];	%high gain edge column vector
    %disp(' test 2kv0')
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = rffutil3(Tmax3,Thi,Tmax4,'hi',angDist);
    
    if strcmp(form,'hf')		%get gain right for pole-zero cancellations
        Tlo = real(Thi) + 1j*imag(c2n(eps));
    else
        Tlo = real(Thi) + 1j*imag(c2n(eps/(w*w)));
    end
    
    
    % The template lies in the  third   quadrant, with the origin on one vertex
    % _________________________________
elseif ( ( (w-wmax>=0) && ((w-wmax)<eps) ) && ( (zmax<= 0) && (zmax > -eps) ) )
    Thi = [ T4(:) ; T1(:) ];	%high gain edge column vector
    %disp(' test 3kv0')
    %change +180 deg to -180 deg to get contiguous phases
    ix = find(real(Thi) > 0);
    if (~isempty(ix)), Thi(ix)= Thi(ix) - 360; end
    Tmax4  = Tmax4 - 360;
    
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = rffutil3(Tmax4,Thi,Tmax1,'hi',angDist);
    
    if strcmp(form,'hf')		%get gain right for pole-zero cancellations
        Tlo = real(Thi) + 1j*imag(c2n(eps));
    else
        Tlo = real(Thi) + 1j*imag(c2n(eps/(w*w)));
    end
    
    
    
    % The template lies in the  fourth   quadrant, with the origin on one vertex
    % _________________________________
elseif ( ( (wmin-w>=0) && ((wmin-w)<eps) )&& ( (zmax<= 0) && (zmax > -eps) ) )
    Thi = [ T1(:) ; T2(:) ];	%high gain edge column vector
    %disp(' test 4kv0')
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = rffutil3(Tmax1,Thi,Tmax2,'hi',angDist);
    
    if strcmp(form,'hf')		%get gain right for pole-zero cancellations
        Tlo = real(Thi) + 1j*imag(c2n(eps));
    else
        Tlo = real(Thi) + 1j*imag(c2n(eps/(w*w)));
    end
    
    
    
    % The origin does not belong to the template
    % _________________________________
    
    
    % The template lies in the first and second quadrants, with the origin outside
    % _________________________________
elseif ( (sign((w-wmin)*(wmax-w)) >= 0) && ( (zmin > eps) ) )
    Thi = [T2(:) ; T3(:) ; T4(:) ];	%high gain edge column vector
    %disp(' test 12kv ')
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = qrff.rffutil3(Tmax2,Thi,Tmax4,'hi',angDist);
    
    Tlo = [T1(:) ];			%low gain edge column vector
    
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Tlo = qrff.rffutil3(Tmax1,Tlo,Tmax1,'lo',angDist);
    
    
    % The template lies in the third and fourth quadrants, with the origin outside
    % _________________________________
elseif ( (sign((w-wmin)*(wmax-w)) >= 0) && ( (zmax <  -eps) ) )
    Thi = [T4(:) ; T1(:) ; T2(:) ];	%high gain edge column vector
    %disp(' test 34kv ')
    %change +180 deg to -180 deg to get contiguous phases
    ix = find(abs(real(Thi)-180) < 1e-10);
    if (~isempty(ix)), Thi(ix)= -180 + 1j*imag(Thi(ix)); end
    ix = find(real(Tmax4)>0);
    if  (~isempty(ix)), Tmax4(ix)=Tmax4(ix)-360; end
    
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = qrff.rffutil3(Tmax4,Thi,Tmax2,'hi',angDist);
    
    Tlo = [T3(:) ];			%low gain edge column vector
    
    % change positive phases to negative
    ix = find(real(Tlo)>0);
    if (~isempty(ix)), Tlo(ix)= Tlo(ix) - 360; end
    ix = find(real(Tmax3)>0);
    if (~isempty(ix)), Tmax3(ix)= Tmax3(ix) - 360; end
    
    Tlo  = qrff.rffutil3(Tmax3,Tlo,Tmax3,'lo',angDist);
    
    
    
    % The template lies in the fourth and first quadrants, with the origin outside
    % _________________________________
elseif ( ( (wmin-w) > eps ) && (sign(zmin*zmax)<=0) )
    Thi = [T1(:) ; T2(:) ; T3(:) ];	%high gain edge column vector
    %disp(' test 41kv ')
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = qrff.rffutil3(Tmax1,Thi,Tmax3,'hi',angDist);
    
    Tlo = [T4(:) ];			%low gain edge column vector
    
    Tlo  = qrff.rffutil3(Tmax4,Tlo,Tmax4,'lo',angDist);	%sort and eliminate and complement
    
    
    % The template lies in the second and third quadrants, with the origin outside
    % _________________________________
elseif (  (wmax-w < -eps)  && (sign(zmin*zmax)<=0) )
    Thi = [T3(:) ; T4(:) ; T1(:) ];	%high gain edge column vector
    %disp(' test 23kv ')
    %add +360 deg to negative phases to get contiguous phases
    ix = find(real(Thi) < 0);
    if (~isempty(ix)), Thi(ix)= Thi(ix) + 360; end
    Tmax1  = Tmax1 + 360;
    
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = qrff.rffutil3(Tmax3,Thi,Tmax1,'hi',angDist);	%
    
    Tlo = [T2(:) ];			%low gain edge column vector
    %add +360 deg to negative phases to get contiguous phases
    ix = find(real(Tlo) < 0);
    if (~isempty(ix)), Tlo(ix)= Tlo(ix) + 360; end
    ix = find(real(Tmax2) < 0);
    if (~isempty(ix)), Tmax2(ix)= Tmax2(ix) + 360; end
    
    
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Tlo  = qrff.rffutil3(Tmax2,Tlo,Tmax2,'lo',angDist);
    
    
    
    % The template lies in the first   quadrant, with the origin outside
    % _________________________________
elseif ( ((wmin-w)  > eps) &&   (zmin > eps)  )
    Thi = [T2(:) ; T3(:) ];	%high gain edge column vector
    %disp(' test 1kv ')
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = qrff.rffutil3(Tmax2,Thi,Tmax3,'hi',angDist);
    
    Tlo = [T1(:) ; T4(:) ];			%low gain edge column vector
    
    Tlo = qrff.rffutil3(Tmax1,Tlo,Tmax4,'lo',angDist);	%sort,  eliminate  doubles,complement
    
    
    
    % The template lies in the second   quadrant, with the origin outside
    % _________________________________
elseif ( ((wmax-w) < -eps) &&   (zmin > eps)  )
    Thi = [T3(:) ; T4(:) ];	%high gain edge column vector
    %disp(' test 2kv ')
    
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = qrff.rffutil3(Tmax3,Thi,Tmax4,'hi',angDist);
    
    Tlo = [T2(:) ; T1(:) ];			%low gain edge column vector
    
    Tlo = qrff.rffutil3(Tmax2,Tlo,Tmax1,'lo',angDist);	%sort etc
    
    
    % The template lies in the third quadrant, with the origin outside
    % _________________________________
elseif ( ((wmax-w) < -eps) &&   (zmax < -eps)  )
    Thi = [T4(:) ; T1(:) ];	%high gain edge column vector
    %disp(' test 3kv ')
    % -360 deg to positive phases to get contiguous phases
    ix = find(real(Thi) > 0);
    if (~isempty(ix)), Thi(ix)= Thi(ix) - 360; end
    ix = find(real(Tmax4)>0);
    if  (~isempty(ix)), Tmax4(ix)=Tmax4(ix)-360; end
    ix = find(real(Tmax1)>0);
    if  (~isempty(ix)), Tmax1(ix)=Tmax1(ix)-360; end
    
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = qrff.rffutil3(Tmax4,Thi,Tmax1,'hi',angDist);
    
    Tlo = [T3(:) ; T2(:) ];			%low gain edge column vector
    
    % -360 deg to positive phases to get contiguous phases
    ix = find(real(Tlo) > 0);
    if (~isempty(ix)), Tlo(ix)= Tlo(ix) - 360; end
    ix = find(real(Tmax3)>0);
    if  (~isempty(ix)), Tmax3(ix)=Tmax3(ix)-360; end
    ix = find(real(Tmax2)>0);
    if  (~isempty(ix)), Tmax2(ix)=Tmax2(ix)-360; end
    
    Tlo = rffutil3(Tmax3,Tlo,Tmax2,'lo',angDist);	%sort etc
    
    
    
    % The template lies in the fourth   quadrant, with the origin outside
    % _________________________________
elseif ( ((wmin-w)  > eps) &&  (zmax < -eps)   )
    Thi = [T1(:) ; T2(:) ];	%high gain edge column vector
    %disp(' test 4kv ')
    %sort, eliminate angular doubles, complement extreme points with phase rounding
    Thi = qrff.rffutil3(Tmax1,Thi,Tmax2,'hi',angDist);
    
    Tlo = [T4(:) ; T3(:) ];			%low gain edge column vector
    
    Tlo = qrff.rffutil3(Tmax4,Tlo,Tmax3,'lo',angDist);	%sort etc
    
    
    
    
    %_______________________________________________________________
end
%---------------------------------------------------------------


ix = find(max( (imag([Tlo,Thi]))' ) < -180);  %eliminate essentially zero

if (~isempty(ix))
    Tlo(ix)=[];
    if isempty(Tlo), Tlo = 0 -1j*313; end
    Thi(ix) = [];
    if isempty(Thi), Thi = 0 -1j*313; end
end

T = [Tlo; Thi];			%upper half=low gain edge,

if strcmp(pzf,'p')		%pole
    T= - [Thi Tlo];
    [~,ix]=sort(real(T(:,1)));	%sort phases in ascending order
    T  = T(ix,:);
    T = T(:);
end

%keyboard




