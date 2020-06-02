function [obj] = rsrs( Tr,M,Ts,Td,w,wco,ordr,Ks,tf,plt,dt,n )
%RSRS    Reference Step Response Specification Calculation
%
%[spc]=rsrs(Tr,M,Ts,Td,w,wco,ordr,Ks,tf,plt,dt,n)
%
%        Transfers time-domain reference step specifications
%        to the frequency domain, using simulations of 2:nd or 3:rd
%        order systems.
%
%        Input variables:
%
%        Tr=[maximum-rise-time, minimum-rise-time, level in percent];
%        Default level is 90%.
%
%        M=overshoot (percent). The default value is M=10%.
%
%        Ts=[Settling time,upper deviation in percent,lower deviation in percent].
%        The default value for the lower deviation is the upper
%        deviation and the default value for the upper deviation is 5%. The
%        default value for the Settling time is 5*Rise-time.
%
%        Td=[Delay-time,level in percent]. Default value for the level is
%        10% and the default value for the delay time is the upper
%        rise-time.
%
%        w=frequency vector. The default value is
%        'logspace(log10(2*pi/(Final-time)/20),log10(2*pi/(Lower
%        rise-time)*20))'.
%
%        wco=cut-off frequency. The lower bound of the frequency
%        specification takes the value '20*log(eps)' for all angular
%        frequencies above wco. The default value is wco=Inf.
%
%        ordr=system order. Can be 2, 3 or 3.1. The default value is 2.
%           order 2 simulates 1/(s^2/w0^2+2*z*s/w0+1).
%           order 3 simulates (s+b)/[(s+a)*(s^2/w0^2+2*z*s/w0+1)]
%                 with appropriate limits on a,b,w0 and z
%                 to avoid resonance peaks, it includes also second order systems
%           order 3.1 is an alternative gridding of order 3, described in
%                 Horowitz: QFT vol 1,  1993. Use instead of 3 for higher speed,
%                   it does not include all the cases of order 2.
%                 See spc_rs2,spc_rs3,spc_rs31 for more information.
%
%        Ks=Static gain. The default value is Ks=1.
%
%        tf=final time. The default value is 2*Setting-time.
%
%        plt=plot option. Logic variable for switching on and off plotting
%        information of step responses. Can be 0 or 1. Default is 1.
%
%        dt=time increment, the default value is 'Final-time/100'.
%
%        n=number of grid-points. The default value is n=[40,40] for the
%        2:nd order system, n=[10,10,10,10] for the 3:rd order
%        system and n=[10,10,10,3] for the alternative gridding. The number
%        of grid-points becomes 'prod(n)'.
%        see the mfiles spc_rs2, spc_rs3, spc_rs31, too see what the different
%        parameterizations are, in order to choose another gridding.
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
%
%
%   Example:
%
%       spc = rsrs([1.2 0.2],10,1.5,[],logspace(-1,2),2.85,3.1);
%
%       Creates a specification SPC, based on simulated reference step 
%       responses that have a rise time of between 0.2 1.2 seconds to the 
%       default level 50%, a maximum overshoot of 10%, a settling time of 
%       1.5 seconds to within 5%. (% refer to % of disturbance step 
%       amplitude). The delay time is the defualt value. The frequency 
%       domain specification is defined for the frequencies logspace(-1,2).
%       A cut-off frequency of 2.85 rad/s is required. The special 
%       simulation system order of 3 is used.
%
%   See also: RSRS, ODSRS, IDSRS, SPC_RS2, SPC_RS3, SPC_RS31.
%

% Author: M Nordin
% Version Upgrade: A. & Y. Greenhut
% OO Version: Daniel Rubin, 7-Jan-2019

if nargin==0
    disp(' qspc = rsrs ( Tr,M,Ts,Td,w,wco,ordr,Ks,tf,plt,dt,n )');
    return
end
if ~(exist('w')==1);  w=[];  end
if ~(exist('dt')==1); dt=[]; end
if ~(exist('n')==1);  n=[];  end
if ~(exist('wco')==1);wco=[]; end
if ~(exist('ordr')==1);ordr=[]; end
if ~(exist('M')==1);M=[]; end
if ~(exist('Ts')==1);Ts=[]; end
if ~(exist('Td')==1);Td=[]; end
if ~(exist('tf')==1);tf=[]; end
if ~(exist('Ks')==1);Ks=[]; end
if ~(exist('plt')==1);plt=[];end
if (length(Tr)==1) || (isempty(Tr))
    error('Both upper and lower rise-times must be given!')
elseif length(Tr)==2
    tr=Tr(1);trl=Tr(2);ptr=90;
elseif length(Tr)==3
    tr=Tr(1);trl=Tr(2);ptr=Tr(3);
end
if trl>tr
    error('Upper Rise-time>=Lower Rise-time must be fulfilled.')
end
if isempty(ordr), ordr=2; end 
if ~isnumeric(ordr)
    error('system order must be a numeric value of either: 2 | 3 | 3.1');
end

if isempty(wco), wco=Inf; end
if isempty(M), M=10; end

if isempty(Ts)
    ts=5*tr; ptsu=5; ptsl=5;
elseif length(Ts)==1
    ts=Ts(1); ptsu=5; ptsl=5;
elseif length(Ts)==2
    ts=Ts(1); ptsu=Ts(2); ptsl=ptsu;
elseif length(Ts)==3
    ts=Ts(1); ptsu=Ts(2); ptsl=Ts(3);
end
if isempty(Td)
    td=tr; ptd=ptr;
elseif length(Td)==1
    td=Td(1); ptd=10;
elseif length(Td)==2
    td=Td(1); ptd=Td(2);
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
    w=logspace(log10(2*pi/tf/20),log10(20*2*pi/trl));
end
e=1000*eps;
tu=[trl,ts,tf].';
stu=[ptr/100,(1+M/100),(1+ptsu/100)].';
tl=[td,tr,ts,tf].';
stl=[0,ptd/100,ptr/100,(1-ptsl/100)].';
if max(~(tu==sort(tu))) || max(~(tl==sort(tl)))
    error('Specification times must fulfill: Delaytime<=Upper Rise-time<=Settingtime and Lower Rise-time<=Setting-time')
end
Iu=diff(tu)>0;
tu=tu(logical([Iu;1]));stu=stu(logical([1;Iu])); %binary sort trick (changed by ilan selig & adi newboer 19/07/2000)
Il=diff(tl)>0;
tl=tl(logical([Il;1]));stl=stl(logical([1;Il]));
tu=sort([0;tu;tu(1:length(tu)-1)+e]);
tl=sort([0;tl;tl(1:length(tl)-1)+e]);
stu=kron(stu,[1;1]);
stl=kron(stl,[1;1]);
t=sort([tu;tl]);
I=[1;diff(t)~=0];
t=t(logical(I));

% Funkar ej !!! tu resp tl måste vara monotont avtagande
stu=interp1(tu,stu,t);
stl=interp1(tl,stl,t);
tab=[t,stu,stl];
switch ordr
    case 2
        [spec_w,spec_t]=qspc.spc_rs2(tab,w,dt,plt,n);
    case 3
        [spec_w,spec_t]=qspc.spc_rs3(tab,w,dt,plt,n);
    case 3.1
        [spec_w,spec_t]=qspc.spc_rs31(tab,w,dt,plt,n);
    otherwise
        error('Invalid system order. Avilable options: 2 (def) | 3 | 3.1');
end
spec_w(:,2:3)=spec_w(:,2:3)+20*log10(Ks);
spec_t(:,2:3)=spec_t(:,2:3)*Ks;
spec_w(spec_w(:,1)>wco,3)=20*log10(eps)*ones(size(spec_w(spec_w(:,1)>wco,3)));

obj = qspc('rsrs',spec_w(:,1),spec_w(:,2),spec_w(:,3),tab,spec_t);

end