function [varargout]=interpnom(obj,w,n)
%INTERPNOM computes nominal plant frequency response by interpolation
%
% Usage: 
%
%   INTERPNOM(P)	computes the nominal case for qplant P using existing
%   templates on a default frequency vector
%
%   INTERPNOM(P,w)     specify the frequency vector 
%
%   INTERPNOM(P,w,n)     also spacifiy the case to be used for interpolation
%           
%   [nom,w]=INTERPNOM(P,...)    returns the nominal case and the frequency 
%   vector, but does not modify P (P.nom is unchangd)
%
% Inputs:
%
%   P           qlant object with precomputed templates
%   w           frequency vector for the intepolated nominal case (defualt
%               is a logarithmic grid of 200 points)
%   n           index of the template case in the 1st template to be used 
%               for the interpolation (default=1)
%
% Remark: if n>1, than the template point #1 is replaced with the selected n
% and all points betewwn 1 and n are shifted by 1 index.

if nargin==0, disp('[nom,w]=INTERPNOM(P,w,n)'), return; end
if nargin<3, n=1; end
if nargin<2, w=[]; end

if isempty(obj.templates)
    error('qplant object must have a nonempty templates property')
end

w_tpl = [obj.templates.frequency];

if isempty(w)
    w=linspace(w_tpl(1),w_tpl(end),200);
    %w=unique([w w_tpl]);
end

tpl=zeros(length(w_tpl),1);
for k=1:length(w_tpl)  
    tpl(k) = obj.templates(k).template(n);
    if n>1
        obj.templates(k).template = [obj.templates(k).template(n) 
            obj.templates(k).template(1:n-1) 
            obj.templates(k).template(n+1:end) ];
        obj.templates(k).parameters = [obj.templates(k).parameters(n),...
            obj.templates(k).parameters(1:n-1),... 
            obj.templates(k).parameters(n+1:end) ];
    end
end

phase = real(tpl);
phase = unwrap(phase*pi/180)*180/pi;
tpl = phase + 1j*imag(tpl);

T = interp1(w_tpl,tpl,w,'linear','extrap'); % 

if nargout==0
    obj.nominal = qfr(T,w);
elseif nargout > 0 
    varargout{1} = T;
end

if nargout>1
    varargout{2}=w; 
end

end

