function [spec_w,spec_t,tab]=spc_od2(spc_tab,w,dt,plt,n)
%SPC_OD2    specification calculation for 2:nd order output disturbance step
% 
%               [spec_w,spec_t,tab]=spc_od2(spc_tab,w,dt,plt,n)
% 
%        Sub-function called by odsrs. Can be used separately for advanced
%        use. The function calculates frequency domain specifications given
%        time domain specifications for a output disturbance reference
%        step. The calculations 
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
%        1-1/(s^2/w0^2+2*z*s/w0+1)
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
%        See also: RSRS, ODSRS, IDSRS, SPC_OD3, SPC_OD31.
% 

% Author: M Nordin
% Version Upgrade: A. & Y. Greenhut
% OO Version: Daniel Rubin, 4-Feb-2019


if nargin==0
    disp(' [spec_w,spec_t,tab] = spc_od2(spc_tab,w,dt,plt,n)')
    return
end
if ~(exist('plt')==1)
    plt=1;
end
if ~(exist('n')==1)
    n=[];
end
if isempty(n)
    n=[40,40];
end
h1=[];h2=[];
ts=spc_tab(:,1);w=w(:);
t=sort([ts(1):dt:ts(length(ts)),ts.']).';
%tab=table1(spc_tab,t); This line has been replaces with the following two
%lines:
spc_tab_size = size(spc_tab); %Y. Greenhut
tab=interp1(spc_tab(:,1),spc_tab(:,2:spc_tab_size(2)),t);

stu=tab(:,1);stl=tab(:,2);
xmin=min(t(stl<stl(1)));
xmax=max(ts);
phimin=eps;
if (min(stl)>-1) && (min(stl)<0)
    phimax=-atan(pi/log(-min(stl)));
elseif min(stl)<-1
    phimax=pi/2-eps;
else
    phimax=eps;
end
qmin=[phimin,xmin];
qmax=[phimax,xmax];
Q=qgrid(n,qmin,qmax);
disp(['Creating array of size 2x',int2str(length(Q(1,:)))]);
z=cos(Q(1,:));
zc=sin(Q(1,:));
%w0=10.^Q(2,:);
w0=2*pi./Q(2,:);

if plt
    figure('Name','Specification calculation','NumberTitle','off');
    hs1=subplot(2,1,1); 
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

T=t;STU=stu;STL=stl;
for j=floor(log2(length(T))):(-1):0
    lz=length(z);
    for k=1:(2^j):length(T)
        st=exp(-w0.*z*T(k)).*(cos(w0.*zc*T(k))+...
            z./zc.*sin(w0.*zc*T(k)));
        I=(st<=STU(k))&(st>=STL(k));
        if plt
            set(h1,'Xdata',T(k)*ones(size(st(I))),'Ydata',st(I),'color','g');
        end
        z=z(I);
        zc=zc(I);
        w0=w0(I);
    end
  if lz~=length(z)
      lz=length(z);
      disp(['Reducing to 2x',int2str(lz)])
  end
  T(1:(2^j):length(T))=[];
  STU(1:(2^j):length(T))=[];
  STL(1:(2^j):length(T))=[];
end

disp(['Number of good disturbance step-responses: ',int2str(length(z))])
if plt
  disp('Plotting...')
end

specu_t=-Inf*ones(size(t));
specl_t=-specu_t;
specu_w=-Inf*ones(size(w));
specl_w=-specu_w;
length(z);
for k=1:length(z)
    st=exp(-w0(k)*z(k)*t).*(cos(w0(k)*zc(k)*t)+...
        z(k)/zc(k)*sin(w0(k)*zc(k)*t));
    mag=20*log10(abs(1-1./(1-w.^2/w0(k)^2+1i*2*z(k)*w/w0(k))));
    specu_t=max(specu_t,st);
    specl_t=min(specl_t,st);
    specu_w=max(specu_w,mag);
    specl_w=min(specl_w,mag);
    if plt
        plot(hs1,t,st);
        semilogx(hs2,w,mag)
    end
end
if plt
    plot(hs1,t,specu_t,'r','linewidth',1.5)
    plot(hs1,t,specl_t,'r','linewidth',1.5)
    semilogx(hs2,w,specu_w,'r','linewidth',1.5)
    semilogx(hs2,w,specl_w,'r','linewidth',1.5)
end
spec_t=[t,specu_t,specl_t];
spec_w=[w,specu_w,specl_w];

end

