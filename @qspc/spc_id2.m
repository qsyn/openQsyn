function [spec_w,spec_t,tab]=spc_id2(spc_tab,w,dt,plt,zmin,n)
%SPC_ID2    Subroutine used by IDSRS 
% 
%               [spec_w,spec_t,tab]=spc_id2(spc_tab,w,dt,plt,zmin,n)
% 
%        Sub-function called by idsrs. Can be used separately for advanced
%        use. The function calculates frequency domain specifications given
%        time domain specifications for an input disturbance step. The
%        calculations are done by gridding the parameters of a second order
%        system. 
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
%        dt is the time increment. (If you want to have a denser or sparser time gridding
% 
%        plt is a logic variable for switching on and off plotting. The
%        default value is 1.
% 
%        zmin is the minimum relative damping of the system. The default
%        value is 1/sqrt(2).        
% 
%        n is the number of grid-points. The default value is n=[10,10,40]
%        and the number of grid-points becomes 'prod(n)'. The first grid
%        variable corresponds to 'phi', the second, to 'x', and the
%        third to 'log10(c)' where the
%        second order system is given by
% 
%        c*s/(s^2+2*z*w0*s+w0^2)
% 
%        where w0=2*pi/x, and z=cos(phi).
%        This is the transfer function corresponding to the input
%        disturbance transfer function when the plant is a general first
%        order system and the controller is a pure integration.
% 
%        Output variables:
% 
%        spec_w; Two-column matrix, where the first column contains the
%        frequency vector and the second column contains the upper bound.
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
%        See also: rsrs, odsrs, idsrs, spc_id2.
% 

% Author: M Nordin
% Version Upgrade: A. & Y. Greenhut
% OO Version: Daniel Rubin, 4-Feb-2019


if nargin==0
  disp('[spec_w,spec_t,tab] = spc_id2(spc_tab,w,dt,plt,zmin,n)')
  return
end
if ~(exist('zmin')==1)
  zmin=1/sqrt(2);
end
if ~(exist('plt')==1)
  plt=1;
end
if ~(exist('n')==1)
  n=[];
end
if isempty(n)
  n=[10,10,30];
end
ts=spc_tab(:,1);w=w(:);
t=sort([ts(1):dt:ts(length(ts)),ts.']).';
spc_tab_size = size(spc_tab);
tab=interp1(spc_tab(:,1),spc_tab(:,2:spc_tab_size(2)),t);
stu=tab(:,1);stl=tab(:,2);
td=min(t(stu==max(stu)));
phimin=sqrt(eps);
phimax=acos(zmin)-sqrt(eps);
xmin=2*pi*td*exp(phimax/tan(phimax));
xmax=max(ts);
ymin=log10(min(stu(length(stu)),-stl(length(stl)))*2*pi/xmax/exp(1));
ymax=log10(2*pi/xmin*exp(1)*max(stu));
qmin=[phimin,xmin,ymin];
qmax=[phimax,xmax,ymax];
Q=qgrid(n,qmin,qmax);
disp(['Creating array of size 3x',int2str(length(Q(1,:)))]);
phi=Q(1,:);
z=cos(Q(1,:));
zc=sin(Q(1,:));
w0=2*pi./Q(2,:);
c=10.^Q(3,:);

if plt
    figure('Name','Specification calculation','NumberTitle','off');
    hs1 = subplot(2,1,1); 
    plot(t,[stu,stl],'r');
    box on
    hold on
    xlabel('time [s]')
    ylabel('Amplitude')
    hs2=subplot(2,1,2);
    set(hs2,'XScale','log')
    hold on
    box on
    axis([w(1) w(length(w)) -30 10]); grid; 
    xlabel('Frequency [rad/s]')
    ylabel('Mag. [dB]')
    drawnow;
end

I=(max(stu)>=c./w0.*exp(-phi.*z./zc))&(min(stl)<=-c./w0.*exp(-(pi+phi).*z./zc))&(c./w0.*exp(-phi.*z./zc)>=min(stu(length(stu)),-stl(length(stl))));
z=z(I);
zc=zc(I);
w0=w0(I);
c=c(I);
disp(['Reducing to 3x',int2str(length(z))])

T=t;STU=stu;STL=stl;
for j=floor(log2(length(T))):(-1):0
  lz=length(z);
  for k=1:(2^j):length(T)
    st=c.*exp(-w0.*z*T(k)).*sin(w0.*zc*T(k))./w0./zc;
    I=(st<=STU(k))&(st>=STL(k));
    c=c(I);
    z=z(I);
    zc=zc(I);
    w0=w0(I);
  end
  if lz~=length(z)
    lz=length(z);
    disp(['Reducing to 3x',int2str(lz)])
  end
  T(1:(2^j):length(T))=[];
  STU(1:(2^j):length(T))=[];
  STL(1:(2^j):length(T))=[];
end

disp(['Number of good step-responses: ',int2str(length(z))])

specu_t=-Inf*ones(size(t));
specl_t=-specu_t;
specu_w=-Inf*ones(size(w));
%mteor=c./w0.*exp(-acos(z).*z./zc);     % unsused!
for k=1:length(z)
  st=c(k).*exp(-w0(k)*z(k)*t).*sin(w0(k)*zc(k)*t)/w0(k)/zc(k);
  mag=20*log10(abs(c(k)*1i*w./(w0(k)^2-w.^2+1i*2*z(k)*w0(k)*w)));
  specu_t=max(specu_t,st); 
  specl_t=min(specl_t,st);
  specu_w=max(specu_w,mag);    
  if plt
    plot(hs1,t,st);
    semilogx(hs2,w,mag)
  end
end
if plt
    plot(hs1,t,specu_t,'r','linewidth',1.5)
    plot(hs1,t,specl_t,'r','linewidth',1.5)
    semilogx(hs2,w,specu_w,'r','linewidth',1.5)
    ylim(ceil([-2 2]+[min(specu_w) max(specu_w)]))
end
spec_t=[t,specu_t,specl_t];
spec_w=[w,specu_w];

                      
