function [Tmax]=frsrs(tpl_nom,tpl,GP,spec,par_nom,par)

%FRSRS      criterion fcn for reference step response specification
%       
%           [Tmax]=frsrs(tpl_nom,tpl,GP,spec,par_nom,par)
%
%           Subroutine to CBND, BNDUPD
%           Criterion fcn 
%               for 20*log10(max|GP/(1+GP)|/min|GP/(1+GP)|)<spec(1)-spec(2) [dB]
%
% Output variables:
%
% Tmax =        value of the criterion function for the n instances of the 
%           nominal open loop candidates to be tested that are
%           contained in GP. Tmax is a row vector of length n
%           with real elements. 
%           The Horowitz bound is the locus of those GP-candidates, 
%           for which Tmax = 0.
%
% Input Variables:
%
% ALL INPUT VARIABLES ARE GIVEN BY THE CALLING BOUND COMPUTATION FCN 
% 
% tpl_nom = the nominal plant template point, a scalar in Nichols form
%               [deg + j*dB].
%
% tpl  =    a m*n plant template matrix where each of the n columns contains 
%           the same template, i.e. n identical columns (of length m)  are  
%           stacked side by side. Each element is in Nichols form [deg + j*dB].
%           (This means each of the m rows have equal elements
%           This is done to simplify matrix computation, the
%           input variables tpl and GP have the same dimension.)
%          
% GP =      A m*n matrix where each column is constant, and each row 
%               contains the n different nominal open loop candidates 
%           [deg + j*dB].
%
% spec =    the specification vector   = 
%           [upper spec dB, lower spec dB, frequency rad/s]
%           corresponding to the frequency for which bounds are 
%           currently computed.
%     
% par_nom = [], not in use. 
%  
% par =     [], not in use. 


% Author: M Nordin
% Version Upgrade: A. & Y. Greenhut


L=n2c(GP+tpl-tpl_nom);              % open loop template in complex form
Tmax=20*log10(max(abs(L)./abs(1+L)));
Tmin=20*log10(min(abs(L)./abs(1+L)));
Tmax=Tmax-Tmin-spec(1)+spec(2);

