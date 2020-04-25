function [T,t_index]=prune(t,Tacc,upper)
%PRUNE      Removes interior points from a connected template
%           [T,t_index]=prune(t,Tacc,upper);
%
%
%   Outputs:
%
%   T:      pruned template vector, with each element of the form
%           deg + j*dB
%
%   t_index: tindex is the indices of the points in t that are contained
%           in the final pruned set. Useful for parameter handling, see
%           the specification of the template file
%
%
%   Inputs:
%
%
%   t:      template vector to be pruned, with each element of the
%           form deg + j*dB
%
%   Tacc:   [degree_accuracy , dB_accuracy], given 2-norm accuracy in Nichols form.
%           A larger value gives a smoother appearance.
%
%   upper:  upper==1	only gives the upper border of a template with
%           phase uncertainty larger than 360 degrees. 
%           default:	upper=0.
%
%
%   Remarks:
%           Removes interior points of a value set in Nichols form deg+j*dB.
%           The points must be connected, i.e. the union of all ellipses 
%           (x-x(i))^/Tacc(1)^2+(y-y(i))^/Tacc(2)^2 where t(i)=x(i)+j*y(i)
%           must be connected.
%
%           Connected templates are produced by e.g. adgrid, adedge, or rff.
%           PRUNE works also on a non-connected template (e.g. one
%           created with the grid method) but not very well. A
%           non-connected template can be made connected
%           with the help of the command CLTMP
%
%           The algorithm handles both wrapping over the Riemann surfaces 
%           and phase uncertainty larger than 360 degrees.
% 
%           T=prune(t,Tacc,1) only gives the upper border of a template with
%           phase uncertainty larger than 360 degrees. 
%
%
%   Examples:      
%           redge5=getfrom('ex2_1b.tpl','t_w5'); % gettpl can also be used!
%               % get a template out of ex2_1b.tpl which
%               % was created with the Recursive Edge Grid method
%           Tredge5=prune(redge5,[2 2]);
%           plot(Tredge5,'ro')
%               % prune and plot
%
%           grid5=getfrom('ex2_1c.tpl','t_w5');
%           Tgrid5=prune(grid5,[2 2]);
%           plot(Tgrid5,'bx')
%               % the same for a template generated with the
%           % Grid method!
%
%  See also TPLPRUNE  TPLREDUC

% This file is part of the OpenQsyn toolox, distribted under GNU LGPL3 license
% Author: M. Nordin


if nargin==0
   disp('  [T,t_index]=prune(t,Tacc,upper);')
   return;
end
if nargin<2 || isempty(Tacc), Tacc=[5,5]; end % empty input case added (DR 6-June-2016)
if nargin<3, upper=0; end
relacc=25;
if length(t)<3
	t_index=1:length(t);
	T=t;
	return;
end
%%Tacc=Tacc*4; % The diameter, not the radius, stupid ??
maxr=relacc+1; %add 1 to compensate for roundoff errors
t1=floor((rem(real(t(:))+5*360,360))*(relacc/Tacc(1))+1+j*imag(t(:))*(relacc/Tacc(2))); %scale with accuracy
minimag=min(imag(t1));
t1=t1-j*minimag+j;   %set left-down corner to 1,1
[maximag tmpi]=max(imag(t1));         % Vertical Size, and starting point (uppermost)
maxreal=ceil(360*relacc/Tacc(1));    % Horizonthal size, always 360 degree
T=t1(tmpi);% The template values
t_index=tmpi; % The indices
tlength=length(t1);
Index=sparse(real(t1),imag(t1),ones(1,tlength)+i*(1:tlength),maxreal,maximag);
 %% The position in Index corresponds to the value of the point. The actual valuse in Index
 %% is a complex number where the real part is the number of points in the original template
 %% and the complex part is the sum of its index number in the original template
r=real(T(1));
k=imag(T(1));
   [i1,j1]=find(Index(max(r-maxr,1):min(r+maxr,maxreal),max(k-maxr,1):min(k+maxr,maximag)));
   if r<maxr+1, %% left wrapping!
	[i2,j2]=find(Index((maxreal+r-maxr):maxreal,max(k-maxr,1):min(k+maxr,maximag)));
   else
	i2=[];j2=[];
   end;
   if r+maxr>maxreal %% right wrapping
	[i3,j3]=find(Index(1:(r+maxr-maxreal),max(k-maxr,1):min(k+maxr,maximag)));
   else
	i3=[];j3=[];
   end;     
   Tind=[i1+max(r-maxr,1)-1+j*(j1+max(k-maxr,1)-1);i2+r-maxr-1+j*(j2+max(k-maxr,1)-1);i3+maxreal+j*(j3+max(k-maxr,1)-1)]-T(1);      
   Tind=Tind((abs(Tind) <= maxr) & (abs(Tind)>0) );
   t_angles=rem(angle(Tind)-acos(abs(Tind)/maxr)+8*pi,2*pi);	
   [tmp tmpi]=min(t_angles);  
	if isempty(tmp);T=t(t_index);return;end;  % only one point found


T(2)=Tind(tmpi)+T(1);
T(2)=rem(real(T(2))+maxreal-1,maxreal)+1+j*imag(T(2)); % compensate for wrapping
if real(Index(real(T(2)),imag(T(2))))==1;
	t_index(2)=imag(Index(real(T(2)),imag(T(2))));
else
	tmpi=find(t1==T(2));
	t_index(2)=tmpi(1);
end;
for n=2:(2*length(t1)+5)
   r=real(T(n));
   k=imag(T(n));
   [i1,j1]=find(Index(max(r-maxr,1):min(r+maxr,maxreal),max(k-maxr,1):min(k+maxr,maximag)));
   if r<maxr+1, %% left wrapping!
	[i2,j2]=find(Index((maxreal+r-maxr):maxreal,max(k-maxr,1):min(k+maxr,maximag)));	
   else
	i2=[];j2=[];
   end
   if r+maxr>maxreal %% right wrapping
      [i3,j3]=find(Index(1:(r+maxr-maxreal),max(k-maxr,1):min(k+maxr,maximag)));
      if size(i3,1)==1 %transfer i3 & j3 into column vector if find return a row vector 
         i3=i3';j3=j3'; %changed by adi newboer &ilan selig (23/7/2000)
      end
   else
	i3=[];j3=[];
   end;     
   Tind=[i1+max(r-maxr,1)-1+j*(j1+max(k-maxr,1)-1);i2+r-maxr-1+j*(j2+max(k-maxr,1)-1);i3+maxreal+j*(j3+max(k-maxr,1)-1)]-T(n);      
   Tind=Tind((abs(Tind) <= maxr) & (abs(Tind)>0) );
   d_T=T(n-1)-T(n);
   d_T=rem(real(d_T)+maxreal*1.5,maxreal)-maxreal*0.5+i*imag(d_T);
   t_angles=rem(angle(Tind/d_T)-acos(abs(Tind)/maxr)-acos(abs(d_T)/maxr)+8*pi,2*pi);	
   [tmp tmpi]=min(t_angles);  
   T(n+1)=Tind(tmpi)+T(n);  
   T(n+1)=rem(real(T(n+1))+maxreal-1,maxreal)+1+j*imag(T(n+1)); % compensate for wrapping
	if real(Index(real(T(n+1)),imag(T(n+1))))==1
		t_index(n+1)=imag(Index(real(T(n+1)),imag(T(n+1))));
	else
		tmpi=find(t1==T(n+1));
		t_index(n+1)=tmpi(1);
    end
   if T(n+1)==T(2) & T(n)==T(1)
		T=T(1:n);
		t_index=t_index(1:n);
      break   
   end
end
if min(imag(T))>1 & upper==0,   %% Phase uncertainty > 360
	disp('PRUNE: Phase uncertainty larger than 360, or not connected')
	[Tlow t_index1]=prune(real(t)-i*imag(t),Tacc,1);
	T=t(t_index);
	T=[T(:).' real(Tlow(:).')-i*imag(Tlow(:).')];
	t_index=[t_index(:).' t_index1(:).'];
else
	T=t(t_index);  % Take the original exact values!!
end 

