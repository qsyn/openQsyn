function [varargout]=showbnd(qdes,bndname,varargin)
%SHOWBND    plots bounds for selected frequencies on Nichols chart
%
% Usage:
% 
%   SHOWBND(qdesign,bndname)    plot computed HS bounds stored in qsedign
%   object, under given name
%
%   SHOWBND(qdesign,bndname,h)    plot on axisting figure given by handle h
%
%   SHOWBND(qdesign,bndname,h,w)    specify frequencies to show
%   
%   SHOWBND(qdesign,bndname,h,w,col)    specify colors
%            
%   h = showbnd2(...)    return figure handle
%
% Output:   
%   
%   h       handle to the figure
%
% Inputs:
%
%   qdes        qdesign object
%   bndname     name of the bounds to plot, strings, e.g. 'odsrs', 'rsrs'; 
%               use 'dom' to plot domiant bounds
%   h           handle to the figure, if =[] a new figure is invoked.  
%   w           frequencies for which to plot each bound, [] gives all freqs.
%   col         RGB array for bound colors
%
%
% Exmaples:
%
%   h = SHOWBND(des,'odsrs',[1 3 5 8]); plots odsrs bounds at frequencies
%   SHOWBND(des,'rsrs',[10 20 30]);     1,3,5,8. plots rsrs bounds for 
%                                       frequencies 10,20,30 on same figure
%   
%   SHOWBND(des,'dom',[],[],[1 0 0])    plots dominant bounds at all
%                                       possible frequencies in red.
%
   
% Adapted from original QSYN SHOWBND. 
% Authors: M Nordin, Copyright P-O Gutman 2004
% NEW OO Method: Daniel Rubin, 2-Jan-2019

if nargin==0
    disp('fhandle = showbnd2(qdesign,bndname,phandle,w,colors)');
	return
end

% check first input
if ~isa(qdes,'qdesign') || ~isscalar(qdes)
    error('first argument must be a qdesign scalar object');
end
% check all other inputs
p = inputParser;
addRequired(p,'bndname',@(x) any(validatestring(x,{qdes.bnd.name,'dom'})));
addOptional(p,'fhandle',[],@(x) any([isgraphics(x), isempty(x)]));
addOptional(p,'w',[],@(x) validateattributes(x,{'numeric'},{'nonnegative','real'}));
addOptional(p,'color',[],@(x) validateattributes(x,{'numeric'},{'ncol',3,'nonnegative','real'}));
parse(p,bndname,varargin{:});

fhandle = p.Results.fhandle;
w = p.Results.w;
col = p.Results.color;

%defaults
if isempty(fhandle)
    fhandle = figure('Name',[bndname,' bounds'],'NumberTitle','off');
    wrap_on=0;
    %hngrid; 
else
    figure(fhandle)
    hold on
    wrap_on=1;
    axis([get(gca,'xlim') get(gca,'ylim')]);  %POGutman 2004-12-01
end
hold on

if strcmp(bndname,'dom')
    disp('computing dominante bounds')
    bnd = dombnd(qdes);
else
    names = {qdes.bnd.name};
    bnd = qdes.bnd(strcmp(names,bndname));
    if isempty(bnd), error('BNDNAME was not found'); end
end

if isempty(w), w=bnd.w; end    
[~,I] = ismember(w,bnd.w);

%col_array = ['m','c','r','g','b','y']'; %def. color changing array 
%col_array = distinguishable_colors(length(w)); 
col_array = qdes.col(I,:);
%%% plot options transfered as a structure
plotstyle=struct('fill',0,'marker','.','color',col_array,'width',1.5); % default settings
if ~isempty(col), plotstyle.color=col; end

ncol=size(plotstyle.color,1);
if length(w)>ncol
    plotstyle.color=repmat(plotstyle.color,ceil(length(w)/ncol),1);
end

% Main Loop
for k=1:length(w)
    bound=bnd.c{I(k)};
    if isempty(bound(~isnan(bound))); return; end
    %The bound is not empty
    if wrap_on
        xlim=get(gca,'xlim'); %degree-axis
        xlim1=[floor(xlim(1)/360) ceil(xlim(2)/360)];
        new_bnd=[];
        while (xlim1(2)-xlim1(1))>0    
            new_bnd=[new_bnd NaN bound+360*(xlim1(1)+1)];
            xlim1(1)=xlim1(1)+1;
        end
        bound=new_bnd;
        ylim=get(gca,'ylim'); %dB-axis
        text_pos=bound(real(bound)>xlim(1) & real(bound)<xlim(2) &...
            imag(bound)>ylim(1) & imag(bound)<ylim(2));
    else
        text_pos=bound(~isnan(bound));
    end
    %if rolling_on;
    %    color=[col_array(rem(k,6)+1),line_style];
    %end;
    h1=plot(bound,'color',plotstyle.color(k,:),'linewidth',plotstyle.width);
    if ~isempty(text_pos)
        text_pos=text_pos(floor((1+length(text_pos))/2));  % peo moved this 960308
        text(real(text_pos),imag(text_pos),num2str(w(k),'%3.3g'),...
            'color',get(h1,'color'));
    else
        %	disp(['The ',bnd,'  bound for ',num2str(w(l)),...
        %	        ' rad/s is outside current axis']); % peo
    end
    
    %ngrid
    xlabel('Phase [deg]')
    ylabel('Magnitude [db]')
    axis tight
    box on
    if nargout>0
        varargout{1} = fhandle;
    end
    
end

