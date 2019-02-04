function [spec_w,spec_t,tab]=spc_od3(spc_tab,w,dt,plt,n)

%SPC_OD3    specification calculation for 3:rd order output disturbance step
% 
%           [spec_w,spec_t,tab]=spc_od3(spc_tab,w,dt,plt,n)
% 
%        Sub-function called by odsrs. Can be used separately for advanced
%        use. The function calculates frequency domain specifications given
%        time domain specifications for an output disturbance reference
%        step. The calculations
%        are done by gridding the parameters of a third order system.
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
%        n is the number of grid-points. The default value is
%        n=[10,10,10,10] and the number of grid-points becomes
%        'prod(n)'. The first grid variable corresponds to '2*pi/p', the
%        second to '2*pi/b', the third to 'log10(u)' and the fourth to
%        'log10(v)', where the third order step-response is given by
% 
%        alpha*exp(-b*t)+beta*exp(-p*t)*sin(ws*t+Phi)
% 
%        where alpha=(u-v)/sqrt(2)+1, beta=(u+v)/sqrt(2). ws and Phi are
%        the corresponding frequency and phase respectively.
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
%        See also: RSRS, ODSRS, IDSRS, SPC_OD2, SPC_OD31.
% 


% Author: M Nordin
% Copyright P-O Gutman 1996
%
% Version Upgrade: A. & Y. Greenhut

if nargin==0
  disp('[spec_w,spec_t,tab]=spc_od3(spc_tab,w,dt,plt,n)')
  return
end;
if ~(exist('plt')==1),
  plt=1;
end;
if ~(exist('n')==1),
  n=[];
end;
if isempty(n)
  n=[10,10,10,10];
end;
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
ymin=min(t(stl<stl(1)));
ymax=max(ts);
qmin=[xmin,ymin,-1.5,-1.5];
qmax=[xmax,ymax,1,1];
Q=qgrid(n,qmin,qmax);
disp(['Creating array of size 4x',int2str(length(Q(1,:)))]);
p=2*pi./Q(1,:);
b=2*pi./Q(2,:);
u=10.^Q(3,:);
v=10.^Q(4,:);
alpha=1+(u-v)/sqrt(2);
beta=(u+v)/sqrt(2);
a=(-beta.^2.*p.^2.*b+2*alpha.*(alpha-1).*b.^2.*p-b.^3.*alpha.^2)./...
  (((alpha-1).*p.^2+b.^2.*alpha-2*b.*alpha.*p).*beta.^2+alpha.*...
  (alpha-1).*b.^2);
w0=sqrt(beta.^2.*p.^2+2*alpha.*(1-alpha).*b.*p+b.^2.*alpha.^2)./...
  sqrt(beta.^2-(alpha-1).^2);
z=p./w0;
ws=w0.*sqrt(1-z.^2);
xi=b./w0;
eta=b./a;
C=(xi.^2-2*z.*xi+1);
A=(xi.^2+eta-2*z.*xi);
B=(xi+z.*xi.^2+z.*eta-xi.*eta-2*z.^2.*xi)./sqrt(1-z.^2);

I=(z>=sqrt(1/(1+(pi/log(max(-min(min(stl),1-eps),eps))^2))))&(b>w0/3)&(a>0);
z=z(I);
a=a(I);
w0=w0(I);
xi=xi(I);
eta=eta(I);
A=A(I);
B=B(I);
C=C(I);
alpha=alpha(I);
beta=beta(I);
p=p(I);
b=b(I);
u=u(I);
v=v(I);
ws=ws(I);

if plt
  hfigbound=figure('Name',['Specification calculation'],'NumberTitle','off');
  subplot(211);
  h1=plot(t,[stu,stl],'r','erasemode','none');
  subplot(212);
  h2=semilogx(w,zeros(size(w)),'r','erasemode','none');
  axis([w(1) w(length(w)) -30 10]);grid;drawnow
end;

disp(['Reducing to 4x',int2str(length(a))])

T=t;STU=stu;STL=stl;
for j=floor(log2(length(T))):(-1):0
  la=length(a);
  for k=1:(2^j):length(T)
    st=((1-eta).*exp(-b*T(k))+exp(-p*T(k)).*(A.*cos(ws*T(k))+...
      B.*sin(ws*T(k))))./C;
    I=(st<=STU(k))&(st>=STL(k));
    if plt
      set(h1,'Xdata',T(k)*ones(size(st(I))),'Ydata',st(I),'color','g');
    end;
    z=z(I);
    a=a(I);
    w0=w0(I);
    xi=xi(I);
    eta=eta(I);
    A=A(I);
    B=B(I);
    C=C(I);
    alpha=alpha(I);
    beta=beta(I);
    p=p(I);
    b=b(I);
    u=u(I);
    v=v(I);
    ws=ws(I);
  end;
  if la~=length(a)
    la=length(a);
    disp(['Reducing to 4x',int2str(la)])
  end;
  T(1:(2^j):length(T))=[];
  STU(1:(2^j):length(T))=[];
  STL(1:(2^j):length(T))=[];
end;

disp(['Number of good step-responses: ',int2str(length(a))])
disp('Plotting...')

specu_t=-Inf*ones(size(t));
specl_t=-specu_t;
specu_w=-Inf*ones(size(w));
specl_w=-specu_w;
for k=1:length(a)
  st=((1-eta(k)).*exp(-b(k)*t)+exp(-p(k)*t).*(A(k).*cos(ws(k)*t)+...
    B(k)*sin(ws(k)*t)))/C(k);
  mag=20*log10(abs(1-(i*w/a(k)+1)./(i*w/b(k)+1)./(1-w.^2/w0(k)^2+...
    i*2*z(k)*w/w0(k))));
  specu_t=max(specu_t,st);
  specl_t=min(specl_t,st);
  specu_w=max(specu_w,mag);
  specl_w=min(specl_w,mag);
  if plt
    set(h1,'Xdata',t,'Ydata',st,'color','b');
    set(h2,'Ydata',mag,'color','b');
  end;
end;
disp('...done.')
if plt;
	set(h1,'Xdata',t,'Ydata',specu_t,'color','r','linewidth',2);
	set(h1,'Ydata',specl_t);
   set(h2,'Ydata',specu_w,'color','r','linewidth',2);
   set(h2,'Ydata',specl_w);
end;
spec_t=[t,specu_t,specl_t];
spec_w=[w,specu_w,specl_w];

