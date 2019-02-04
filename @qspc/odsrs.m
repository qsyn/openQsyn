function obj = odsrs(Tr,M,Ts,Td,w,ordr,Ks,tf,plt,dt,n)
%ODSRS      computation of output disturbance step response specification
% 
%   obj = ODSRS(Tr,M,Ts,Td,w,ordr,Ks,tf,plt,dt,n)
% 
%        Transfers time-domain output disturbance step specifications
%        to the frequency domain, using simulations of a 2:nd or 3:rd
%        order system and gridding over its parameters.
%
%
%Output
% 
%        spec_w; two-column matrix, where the first column contains the
%        frequency vector (rad/s), the second column contains the upper bound (dB).
% 
%        spec_t; two-column matrix, where the first column contains the
%        time vector (sec), the second column contains the upper envelope of the
%        step responses and the third column contains the lower envelope.
% 
%        tab; Three-column matrix, where the first column contains the time
%        vector corresponding to the time-domain specifications. The second
%        column contains the upper bound corresponding to the time-domain
%        specifications and the third column contains the corresponding
%        lower bound.
% 
% 
%Inputs
%
%   Tr      [Maximum rise-time, Minimum rise-time, level (percent)];
%           Default level is 10%.
% 
%   M       undershoot (percent). The default value is M=10%.
% 
%   Ts      [Settling time, upper deviation (percent),lower deviation (percent)]. 
%           The default value for the lower deviation is the upper deviation 
%           and the default value for the upper deviation is 5%. The default
%           value for the Settling time is 5*Rise-time.
% 
%   Td      [Delay-time,level (percent)]. Default value for the level is
%           90% and the default value for the delay time is the upper
%           rise-time.
% 
%  	w       frequency vector. The default value is
%           'logspace(log10(2*pi/(Final-time)/2),log10(2*pi/(Lower-rise-time)*2))'. 
% 
%   ordr    system order. Can be 2, 3 or 3.1. The default value is 2. If
%           ordr is set to 3.1 an alternative gridding for the 3rd order
%           system is performed.
% 
%   Ks      Initial gain. The default value is Ks=1.
%
%   tf      final time. The default value is 2*Settling-time.
% 
%   plt     plot option. Logic variable for switching on and off plotting
%           information of step responses. Can be 0 or 1. Default is 1.
%
%   dt      time increment; The default value is 'Final-time/100'.
% 
%   n       number of grid-points. The default value is n=[40,40] for the
%           2nd order system, n=[10,10,10,10] for the 3rd order system and
%           n=[10,10,10,3]. For the alternative gridding the number of 
%           grid points becomes 'prod(n)'. 
%           See spc_od2, spc_od3, spc_od31 for more about the parametrization
%
%
%Example:
% 
%	spc1 = odsrs([0.3 0.1 50],50,1.5,[],logspace(-1,2),2);
%       
%       constracts a qspc opject called spec1 based on simulated output 
%       disturbance step responses that have a rise time of between 0.1 0.3
%       seconds to the level 50%, a maximum undershoot of 50%, a settling 
%       time of 1.5 seconds to within 5%. (% refer to % of disturbance step 
%       amplitude). The frequency domain specification is defined for the 
%       frequencies logspace(-1,2). The imulation system order of 2 is used. 
%
% 
%See also: qspc.rsrs qspc.idsrs qspc.idsrs, qspc.SPC_OD2, qspc.SPC_OD3, qspc.SPC_OD31.
% 

% Author: M Nordin
% Version Upgrade: A. & Y. Greenhut
% OO Version: Daniel Rubin, 4-Feb-2019


if nargin==0
  disp('  spc = odsrs(Tr,M,Ts,Td,w,ordr,Ks,tf,plt,dt,n)');
  return
end

if ~(exist('w')==1)    
  w=[];
end
if ~(exist('dt')==1)
  dt=[];
end
if ~(exist('n')==1)
  n=[];
end
if ~(exist('wco')==1)
  wco=[];
end
if ~(exist('ordr')==1)
  ordr=[];
end
if ~(exist('M')==1)
  M=[];
end
if ~(exist('Ts')==1)
  Ts=[];
end
if ~(exist('Td')==1)
  Td=[];
end
if ~(exist('tf')==1)
  tf=[];
end
if ~(exist('Ks')==1)
  Ks=[];
end
if ~(exist('plt')==1)
  plt=[];
end
if (length(Tr)==1) || (isempty(Tr))
  error('Both upper and lower rise-times must be given!')
elseif length(Tr)==2
  tr=Tr(1);trl=Tr(2);ptr=10;
elseif length(Tr)==3
  tr=Tr(1);trl=Tr(2);ptr=Tr(3);
end
if trl>tr
  error('Upper Rise-time>=Lower Rise-time must be fulfilled.')
end
if isempty(ordr)
  ordr=2;
end
if ~isnumeric(ordr)
    error('system order must be a numeric value of either: 2 | 3 | 3.1');
end
if isempty(wco)
  wco=0;
end
if isempty(M)
  M=10;
end
if isempty(Ts)
  ts=5*tr;ptsu=5;ptsl=5;
elseif length(Ts)==1
  ts=Ts(1);ptsu=5;ptsl=5;
elseif length(Ts)==2
  ts=Ts(1);ptsu=Ts(2);ptsl=ptsu;
elseif length(Ts)==3
  ts=Ts(1);ptsu=Ts(2);ptsl=Ts(3);
end
if isempty(Td)
  td=tr;ptd=ptr;
elseif length(Td)==1
  td=Td(1);ptd=90;
elseif length(Td)==2
  td=Td(1);ptd=Td(2);
end
if isempty(tf)
  tf=2*ts;
end
if isempty(Ks)
  Ks=1;
end
if isempty(plt)
  plt=1;
end
if isempty(dt)
  dt=tf/100;
end
if isempty(w)
  w=logspace(log10(2*pi/tf/2),log10(2*2*pi/trl));
end
e=1000*eps;
tl=[trl,ts,tf].';
stl=[ptr/100,-M/100,-ptsl/100].';
tu=[td,tr,ts,tf].';
stu=[1,ptd/100,ptr/100,ptsu/100].';
if max(~(tl==sort(tl))) || max(~(tu==sort(tu)))
  error('Specification times must fulfill: Delaytime<=Upper Rise-time<=Setting-time and Lower Rise-time<=Setting-time')
end
Il=diff(tl)>0;
tl=tl(logical([Il;1]));stl=stl(logical([1;Il]));
Iu=diff(tu)>0;
tu=tu(logical([Iu;1]));stu=stu(logical([1;Iu]));
tl=sort([0;tl;tl(1:length(tl)-1)+e]);
tu=sort([0;tu;tu(1:length(tu)-1)+e]);
stl=kron(stl,[1;1]);
stu=kron(stu,[1;1]);
t=sort([tl;tu]);
I=[1;diff(t)~=0];
t=t(logical(I));
stl=interp1(tl,stl,t);
stu=interp1(tu,stu,t);
tab=[t,stu,stl];
switch ordr
    case 2
        [spec_w,spec_t] = qspc.spc_od2(tab,w,dt,plt,n);
    case 3
        [spec_w,spec_t] = spc_od3(tab,w,dt,plt,n);
    case 3.1
        [spec_w,spec_t] = qspc.spc_od31(tab,w,dt,plt,n);
    otherwise
        error('Invalid system order. Avilable options: 2 (def) | 3 | 3.1');
end
spec_w(:,2:3)=spec_w(:,2:3)+20*log10(Ks);
spec_t(:,2:3)=spec_t(:,2:3)*Ks;
spec_w(spec_w(:,1)<wco,3)=20*log10(eps)*ones(size(spec_w(spec_w(:,1)<wco,3)));
%insert(spcfile,spec_w(:,1:2),[specname,'_w'],'r'); 
%insert(spcfile,spec_t,[specname,'_t'],'r'); 
%insert(spcfile,tab,[specname,'_tab'],'r'); 

obj = qspc('odsrs',spec_w(:,1),spec_w(:,2),[],tab,spec_t);

end
