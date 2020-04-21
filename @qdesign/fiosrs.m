function [Tmax]=fiosrs(tpl_nom,tpl,GP,spec,par_nom,par)

%FIOSRS     criterion fcn for input step response spec/complementary sensitivity
%           [Tmax]=fiosrs(tpl_nom,tpl,GP,spec,par_nom,par)
%
%           subroutine for CBND, BNDUPD
%           criterion function for spec(2) < |GP/(1+GP)| < spec(1)
%           criterion function also for 1 d-o-f servo specifications
%
%   Output:
%
%   Tmax        value of the criterion 
%
%
%   Inputs:
%
% ALL INPUT VARIABLES ARE GIVEN BY THE CALLING BOUND COMPUTATION FCN
%
%   tpl_nom     plant template nominal, Pnom, scalar in Nichols form, deg+j*dB  
%
%   tpl         plant template P, in Nichols form
%
%   GP          open loop nominal candidate, Lnom=G*Pnom, in Nichols form, deg+j*dB
%
%   spec        specification values = [upper limit dB, lower limit dB, frequency rad/s], 
%		which equals the specification row and the appropriate frequency, 
%               of the current
%               specification.
%
%   par_nom     [], not in use
%
%   par         [], not in use



% Author: M Nordin
% Version Upgrade: A. & Y. Greenhut

L=n2c(GP+tpl-tpl_nom);
Tmax=20*log10(max(abs(L)./abs(1+L)));
Tmax=max([Tmax-spec(1);spec(2)-Tmax]); 
% must fulfill both upper and lower specification
