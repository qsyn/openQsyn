function [spec_w,spec_t,tab]=spc_rs2(spc_tab,w,dt,plt,n)
%SPC_RS2    specification calculation for 2:nd order reference step
% 
%           [spec_w,spec_t,tab]=spc_rs2(spc_tab,w,dt,plt,n)
% 
%        Sub-function called by rsrs. Can be used separately for advanced
%        use. The function calculates frequency domain specifications given
%        time domain specifications for a reference step. The calculations
%        are done by gridding the parameters of a second order system.
% 
%        Input variables:
% 
%        spc_tab is a 3-column matrix containing time domain specifications
%        for a reference step. The first column contains the times, second
%        column contains the upper bound and the third contains the lower
%        bounds. The upper and lower bounds are assumed to be LINEARLY
%        interpolated between the values in the second and third columns,
%        respectively. 
% 
%        w is the frequency vector, preferably created by logspace.
% 
%        dt is the time increment.
% 
%        plt is a logic variable for switching on and off plotting. The
%        default value is 1.
% 
%        n is the number of grid-points. The default value is n=[40,40]
%        and the number of grid-points becomes 'prod(n)'. The first grid
%        variable corresponds to 'phi' and the second, to 'x', where the
%        second order system is given by
% 
%        1/(s^2/w0^2+2*z*s/w0+1)
% 
%        where w0=10^x, and z=cos(phi).
% 
%        Output variables:
% 
%        spec_w; Three-column matrix, where the first column contains the
%        frequency vector, the second column contains the upper bound and
%        the third column contains the lower bound.
% 
%        spec_t; Three-column matrix, where the first column contains the
%        time vector, the second column contains the upper envelope of the
%        step responses and the third column contains the lower envelope.
% 
%        tab; Three-column matrix, where the first column contains the time
%        vector corresponding to the time-domain specifications. The second
%        column contains the upper bound corresponding to the time-domain
%        specifications and the third column contains the corresponding
%        lower bound.
% 
%        See also: RSRS, ODSRS, IDSRS, SPC_RS3, SPC_RS31.
% 


% Author: M Nordin
% Version Upgrade: A. & Y. Greenhut

if nargin==0
  disp('[spec_w,spec_t,tab]=spc_rs2(spc_tab,w,dt,plt,n)')
  return % instead of break (D.R. 11/7/15)
end;
if ~(exist('plt')==1)
  plt=1;
end;
if ~(exist('n')==1)
  n=[];
end;
if isempty(n)
  n=[40,40];
end;
h1=[];h2=[];
ts=spc_tab(:,1);w=w(:);
t=sort([ts(1):dt:ts(length(ts)),ts.']).';

% Funkar ej !!! spc_tab måste vara monotont avtagande
spc_tab_size = size(spc_tab); %Y. Greenhut
tab=interp1(spc_tab(:,1),spc_tab(:,2:spc_tab_size(2)),t);
stu=tab(:,1);stl=tab(:,2);
xmin=min(t(stu>stu(1)));
xmax=max(ts);
phimin=eps;
if (max(stu)<2)&(max(stu)>1),  
  phimax=-atan(pi/log(max(stu)-1));
elseif max(stu)>2
  phimax=pi/2-eps;
else
  phimax=eps;
end;
qmin=[phimin,xmin];
qmax=[phimax,xmax];
Q=qgrid(n,qmin,qmax);
disp(['Creating array of size 2x',int2str(length(Q(1,:)))]);
z=cos(Q(1,:));
zc=sin(Q(1,:));
w0=2*pi./Q(2,:);
if plt
  hfigbound=figure('Name',['Specification calculation'],'NumberTitle','off');
  subplot(211);
  h1=plot(t,[stu,stl],'r','erasemode','none');
  subplot(212);
  h2=semilogx(w,zeros(size(w)),'r','erasemode','none');
  axis([w(1) w(length(w)) -30 10]);grid;drawnow
end;

T=t;STU=stu;STL=stl;
for j=floor(log2(length(T))):(-1):0
  lz=length(z);
  for k=1:(2^j):length(T)
    st=1-exp(-w0.*z*T(k)).*(cos(w0.*zc*T(k))+...
    z./zc.*sin(w0.*zc*T(k)));   
    I=(st<=STU(k))&(st>=STL(k));
    if plt
      set(h1,'Xdata',T(k)*ones(size(st(I))),'Ydata',st(I),'color','g');
    end;
    z=z(I);
    zc=zc(I);
    w0=w0(I);
  end;
  if lz~=length(z)
    lz=length(z);
    disp(['Reducing to 2x',int2str(lz)])
  end;
  T(1:(2^j):length(T))=[];
  STU(1:(2^j):length(T))=[];
  STL(1:(2^j):length(T))=[];
end;

disp(['Number of good step-responses: ',int2str(length(z))])
if plt
  disp('Plotting...')
end;

specu_t=-Inf*ones(size(t));
specl_t=-specu_t;
specu_w=-Inf*ones(size(w));
specl_w=-specu_w;

for k=1:length(z)
  st=1-exp(-w0(k)*z(k)*t).*(cos(w0(k)*zc(k)*t)+...
  z(k)/zc(k)*sin(w0(k)*zc(k)*t));
  mag=20*log10(abs(1./(1-w.^2/w0(k)^2+i*2*z(k)*w/w0(k))));
  specu_t=max(specu_t,st); 
  specl_t=min(specl_t,st);
  specu_w=max(specu_w,mag);
  specl_w=min(specl_w,mag);
  if plt
    set(h1,'Xdata',t,'Ydata',st,'color','b');
    set(h2,'Ydata',mag,'color','b');
  end;   
end;
if plt
	disp('...done.')
	set(h1,'Xdata',t,'Ydata',specu_t,'color','r','linewidth',2);
	set(h1,'Ydata',specl_t);
   set(h2,'Ydata',specu_w,'color','r','linewidth',2);
   set(h2,'Ydata',specl_w);
end;
spec_t=[t,specu_t,specl_t];
spec_w=[w,specu_w,specl_w];
 

