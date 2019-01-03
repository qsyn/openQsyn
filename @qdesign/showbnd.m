function [varargout]=showbnd(qdes,bndname,phandle,w,colors)
%SHOWBND    plots bounds from a QDESIGN object for selected frequencies, on
%Nichols chart
%
%   Usage:   showbnd2(qdesign,bndname)
%            showbnd2(qdesign,bndname,phandle)
%            showbnd2(qdesign,bndname,phandle,w)
%            showbnd2(qdesign,bndname,phandle,w,colors)
%            phandle = showbnd2(qdesign,bndname,phandle,w,colors)
%   
%   Output:
% 
%   phandle     handle to the figure  
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

if nargin<2, error('Not enough input arguments!'); end
if ~exist('phandle'), phandle=[]; end
if ~exist('w'), w=[]; end
if ~exist('colors'), colors=[]; end

%defaults
if isempty(phandle)
    phandle = figure('Name',[bndname,' bounds'],'NumberTitle','off');
    wrap_on=0;
    %hngrid; % WHY NOT??? (Daniel R 19-June-2016)
else
    figure(phandle)
    wrap_on=1;
    axis([get(gca,'xlim') get(gca,'ylim')]);  %POGutman 2004-12-01
end
hold on;

names = [qdes.bnd.name];
bnd = qdes.bnd(strcmp(names,bndname));

if isempty(w), w=bnd.w; end    

col_array = ['m','c','r','g','b','y']'; %def. color changing array 
%%% plot options transfered as a structure
plotstyle=struct('fill',0,'marker','.','color',col_array,'width',1.5); % default settings
if ~isempty(colors), plotstyle.color=colors; end

ncol=size(plotstyle.color,1);
if length(w)>ncol
    plotstyle.color=repmat(plotstyle.color,ceil(length(w)/ncol),1);
end

% Main Loop
for k=1:length(w)
    bound=bnd.c{k};
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

    if nargout>0
        varargout{1}=h1;
    end
    
end

