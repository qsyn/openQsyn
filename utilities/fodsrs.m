function[Smax]=fodsrs(tpl_nom,tpl,GP,spec,par_nom,par)

%FODSRS  	criterion function for output dist step resp specification or sensitivity
%           function[Smax]=fodsrs(tpl_nom,tpl,GP,spec,par_nom,par)
%
%           Subroutine to CBND, BNDUPD
%           Criterion function for |1/(1+GP)| < spec(1)
%
%   Output:
%
%   Smax        value of the criterion 
%
%
%   Inputs:
%
% ALL INPUT VARIABLES ARE GIVEN BY THE CALLING BOUND COMPUTATION FCN
%
%
%   tpl_nom     plant template nominal, Pnom, scalar in Nichols form, deg+j*dB  
%
%   tpl         plant template P, in Nichols form
%
%   GP          open loop nominal candidate, Lnom=G*Pnom, in Nichols form, deg+j*dB
%
%   spec        [sensitivity specification value dB, frequency rad/s]
%
%   par_nom     [], not in use
%
%   par         [], not in use

% Author: P-O Gutman, M Nordin
% Copyright: El-Op Electro-Optics Industries Ltd, and Novomatic AB 1996
% Version Upgrade: A. & Y. Greenhut

Smax=-20*log10(min(abs(1+n2c(GP+tpl-tpl_nom))))-spec(1); %obs!	

