function [PSmax] = fidsrs(tpl_nom,tpl,GP,spec,par_nom,par)

%FIDSRS     criterion function for input disturbance step response specification
%           [PSmax]=fidsrs(tpl_nom,tpl,GP,spec,par_nom,par)
%
%           subroutine for CBND, BNDUPD
%           criterion function for  |P/(1+GP)| < spec(1) 
%
%   Output:
%
%   PSmax        value of the criterion 
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
%   spec        [specification value dB, frequency rad/s] which equals
%               the appropriate row, and  frequency, of the current
%               specification.
%
%   par_nom     [], not in use
%
%   par         [], not in use
%   
% 

% Author: M Nordin
% Version Upgrade:A. &  Y. Greenhut

L=n2c(GP+tpl-tpl_nom);
PSmax=20*log10(max(abs(n2c(tpl))./abs(1+L)))-spec(1);

