function [varargout]=show(obj ,varargin)
%SHOW    displays a template object on Nichols chart 
%     
%  	show(QTPL)   displays template QTPL 
%   QTPL.show    alternative usage  
%
%   show(TPLf,FHAND)    draws the Nichols chart in figure with handle FHAND   
%
%   show(TPLF,PARAMETER,VALUE)   use parameter/value pairs to
%   specify additional properties:
%       use PARAMETER='color'   with VALUE a color array in RGB format
%       use PARAMETER='marker'  with VALUE a string for marker points 
%       use PARAMETER='fill'    with VALUE equals 1 (filled) or zero (def) 
%       use PARAMETER='case'    with VALUE a vector of indices specifing
%                               plant case(s) to show
%
%   H=showtpl(...)      returns ahandle to the figure 
%
%
%   Remarks:            During the execution of the display under the
%                   'point' option, pressing  the left button gives the next point and 
%                   pressing the right button terminates showtpl
% 
%                   If  phandle is a figure with no axes (e.g. after
%                   you have opened a new figure with the command
%                   figure), then PLEASE issue the command hngrid
%                   before showtpl, so that the axes will be
%                   suitable for the display of the templates.
% 
%                   You may zoom in on details of the template display
%                   by either the Matlab command zoom, or by the
%                   Qsyn command hzoom, that puts a zoom menu on the
%                   current figure toolbar, that also lets you zoomout
%                   beyond the borders of the original picture.
% 
%   Examples:       TODO!!
% 
%               
%   See also         HNGRID MGRID CTPL HZOOM

% Adapted from QSYN SHOWTPL. 
% Authors: B. Cohen, P-O Gutman and M Nordin, 1996
% Version upgrade: A. & Y. Greenhut
% NEW OO Method: Daniel Rubin, 26-Dec-2018

if nargin==0, disp('  h = show(qtpl,...)'), return; end
%if ~isscalar(tpl), error('QTPL.SHOW does not support QTPL arrays. use ')
    
phandle=[];
IDX=[]; % defaults

col_array = ['m','c','r','g','b','y']'; %def. color changing array 
plotstyle=struct('fill',0,'marker','.','color',col_array); % default settings

if nargin>1
    k=1;
    if ishandle(varargin{1}), phandle=varargin{1}; k=2; end
    if isempty(varargin{1}), k=2; end
    
    while k<=nargin-2
        if ischar(varargin{k})
            if nargin<(k+1), error('incorrect parameter/value pairing!'); end
            PARAMETER=varargin{k};
            switch PARAMETER
                case 'color', plotstyle.color=varargin{k+1};
                case 'marker', plotstyle.marker=varargin{k+1};
                case 'fill', plotstyle.fill=varargin{k+1};
                case 'case', IDX=varargin{k+1};
                otherwise, error('unknown parameter %s',PARAMETER);
            end
        else
            error('incorrect parameter/value pairing!');
        end
        k=k+2;
    end
end

N = numel(obj);
w_op = [obj.frequency]';

% Allow numeric 'case' input to specify the plant 
if ~isnumeric(IDX), error('IDX must be a numeric array'); end

ncol=size(plotstyle.color,1);
if length(w_op)>ncol
    plotstyle.color=repmat(plotstyle.color,ceil(length(w_op)/ncol),1);
end
Lstyle=plotstyle.marker;

if isempty(phandle)
    phandle = figure('Name','Templates on Nichols Chart','NumberTitle','off');
end
    
hold on
for k=1:N
    if isempty(IDX)
        tpl = obj(k).template;
    else
        tpl = obj(k).template(IDX);
    end
    ntpl=size(tpl,1);
    Lcolor = plotstyle.color(k,:);
    if plotstyle.fill
        patch(real(tpl),imag(tpl),Lcolor,'erasemode','xor');
    end
    plot(real(tpl),imag(tpl),'Color',Lcolor,'LineStyle','none',...
        'Marker',Lstyle,'Markersize',4,'zdata',1:ntpl) 
    plot(real(tpl(1)),imag(tpl(1)),'s','Color',Lcolor)     
    text(real(tpl(1)+1),imag(tpl(1))-1,[num2str(w_op(k))]);
    
end

% pause(0.01); drawnow 
% XLIM=xlim;
% YLIM=ylim;
% hngrid;
% axis([XLIM YLIM]);

dcm = datacursormode(gcf);
set(dcm,'updatefcn',@datatipfunc)

if nargout==1
    varargout{1}=phandle;
end
 
function output_txt = datatipfunc(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
%dts = get(event_obj.Target,'zdata');
output_txt = {['Phase: ',num2str(pos(1),4), '[deg]'],...
              ['Mag: ',num2str(pos(2),4),' [db]'],...
              ['par set: ',num2str(pos(3))]
              };

