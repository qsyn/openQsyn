function obj = lag(Freq,beta)
%LAG returns a QCTRL lag compensator
%
% Usage:
%
% clag = qctl.lag(Phase,W)  returns a 1st oreder lag compensator
% with a given phase-lag at frequeny W.
%
if nargin==2 % order=1
      s = qctrl(0,[],1);
      obj = (10*s+Freq)/(10*s+Freq/beta);
else
    error('incorrect number of arguments')
end
end