function [varargout]=showbnd(qdes,bndname,fhandle,w,colors)
%SHOWBND    plots bounds from a QDESIGN object for selected frequencies, on
%Nichols chart
%
%   Usage:   showbnd2(qdesign,bndname)
%            showbnd2(qdesign,bndname,phandle)
%            showbnd2(qdesign,bndname,phandle,w)
%            showbnd2(qdesign,bndname,phandle,w,colors)
%            fhandle = showbnd2(qdesign,bndname,phandle,w,colors)
%   
%   Output:
% 
%   fhandle     handle to the figure  
%
%   Inputs:
%
%   bndfile     the boundfile name, string with or without extension '.bnd',
%               containing the bounds
%    
%   phandle     handle to the figure, if =[] a new figure is invoked.  
%
%   bnd         name of the bounds to plot, strings, e.g. 'odsrs', 'rsrs'
%
%   w           frequencies for which to plot each bound, [] gives   all freqs.
%
%   colors      RGB array for bound colors

% Adapted from QSYN SHOWBND. 
% Authors: M Nordin, Copyright P-O Gutman 2004
% NEW OO Method: Daniel Rubin, 2-Jan-2019

if nargin==0
    disp('showbnd(bndfile,bnd,w,phandle,colors)');
	return
end

if nargin<5, colors=[]; end
if nargin<4, w=[]; end
if nargin<3, fhandle=[]; end
if nargin<2, error('Not enough input arguments!'); end

%defaults
if isempty(fhandle)
    fhandle = figure('Name',[bndname,' bounds'],'NumberTitle','off');
    wrap_on=0;
    %hngrid; % WHY NOT??? (Daniel R 19-June-2016)
else
    figure(fhandle)
    hold on
    wrap_on=1;
    axis([get(gca,'xlim') get(gca,'ylim')]);  %POGutman 2004-12-01
end
hold on;

names = {qdes.bnd.name};
bnd = qdes.bnd(strcmp(names,bndname));
if isempty(bnd), error('BNDNAME was not found'); end

if isempty(w), w=bnd.w; end    

%col_array = ['m','c','r','g','b','y']'; %def. color changing array 
col_array = distinguishable_colors(length(w)); 
%%% plot options transfered as a structure
plotstyle=struct('fill',0,'marker','.','color',col_array,'width',1.5); % default settings
if ~isempty(colors), plotstyle.color=colors; end

ncol=size(plotstyle.color,1);
if length(w)>ncol
    plotstyle.color=repmat(plotstyle.color,ceil(length(w)/ncol),1);
end

[~,I] = ismember(w,bnd.w);

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
    
    ngrid
    xlabel('Phase [deg]')
    ylabel('Magnitude [db]')
    axis tight
    if nargout>0
        varargout{1} = fhandle;
    end
    
end

