function obj = lag(Freq,beta)
%LAG returns a QCTRL lag compensator
%
% Usage:
%
% clag = QCTRL.LAG(Phase,W)  returns a 1st oreder lag compensator
% with a given phase-lag at frequeny W.
%
if nargin==2 % order=1
      s = qctrl(0,[],1);
      obj = (10*s+2)/(10*s+2/3);
else
    error('incorrect number of arguments')
end
end