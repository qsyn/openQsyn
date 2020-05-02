function [varargout]=show(obj,plot_op,linespec,phandle)
%SHOWSPC    displays a single specification 
%               
% [phandle_new] = showspc(spcObject,plot_op,linespec,phandle)
%
%   OUTPUTS:
%   
%       phandle_new:    A handle to the figure.
%
%   INPUTS:
%
%       plot_op     'both' (def) plot both frequncy and time domain specs. 
%                   'time' plots only time domain specs.
%                   'freq' plots only frequency domain specs.
% 
%       linespec    line specification, e.g. '-r' (def), 'b--',....
%
%       phandle     A handle to the figure you want to plot in, use 
%                   phandle=gcf if you want to plot in the current fig.
%                   An empty phandle invokes a new figure (default).
% 
%                   You may zoom in on details of the template display
%                   by either the Matlab command zoom, or by the
%                   Qsyn command hzoom, that puts a zoom menu on the
%                   current figure toolbar, and which also lets you zoom out
%                   beyond the borders of the original picture.
%
%   Example:
%               
%   SHOW(spec,'freq','g',gcf)   plots the frequency domain odsrs specs. 
%   from a QSPC object spec in the current figure (gcf) in green.
%
%
%   See also       RSRS, ODSRS, SPCUPD, IDSRS

% Adapted from QSYN SHOWSPC. 
% Authors: M Nordin,  A. & Y. Greenhut,  Copyright P-O Gutman 1996
% NEW OO Method: Daniel Rubin, 7-Jan-2019

if nargin==0
   disp('  [phandle]=show(qspc,plot_op,color,phandle)')
   return
end

if ~exist('phandle','var'), phandle=[]; end
if ~exist('linespec','var'), linespec=[]; end
if ~exist('plot_op','var'), plot_op=[]; end

if isempty(plot_op); plot_op='both'; end
if isempty(linespec); linespec='r-'; end

if isempty(phandle) %new figure
	phandle_new=figure('Name',['specifications for ',obj.name],'NumberTitle','off');
else
	figure(phandle)
    phandle_new = phandle; 
	axis(axis); hold on
end

if isempty(obj.timeres)
    plot_op = 'freq';
else
    spec_t = obj.timeres;
    tab = obj.timespc;
end

if strcmp(plot_op,'both') 
    subplot(2,1,1)
end

if strcmp(plot_op,'time') || strcmp(plot_op,'both')
    plot(spec_t(:,1),spec_t(:,2),linespec,'linewidth',2); hold on;
    plot(spec_t(:,1),spec_t(:,3),linespec,'linewidth',2);
    plot(tab(:,1),tab(:,2:3),linespec);
    xlabel('Time [s]')
    ylabel('normeliized response')
end

if strcmp(plot_op,'both') 
    subplot(2,1,2)
end

if strcmp(plot_op,'freq') || strcmp(plot_op,'both')
    semilogx(obj.frequency,obj.upper,linespec); 
    if ~isempty(obj.lower)
        hold on
        semilogx(obj.frequency,obj.lower,linespec); 
    end
    ylim(round([min(obj.upper)-5 max(obj.upper)+2]))
    xlim([obj.frequency(1) obj.frequency(end)])
    xlabel('Frequency [rad/s]')
    ylabel('Spec. [db]')
end
    
if nargout==1
    varargout={phandle_new};
end

end
    
