function qngrid(varargin)
%QNGRID     draws T=L/(1+L) or S=1/(1+L) loci in a Nichols diagram
%
%   Usage: 
%   
%   QGRID()  draws default  T=L/(1+L) loci
%   
%   QGRID(Type)  specifiy type: T (def) | S
%
%   QGRID(Type,ConstMag,ConstPhase)  specifiy also values for the loci contours 
%
%   QGRID(...,parameter,value)  additional
%   parameter/value pairs
%
% 
%   Optional Inputs:
%
%           Type:       'T' (default)draws L/(1+L) loci | 'S' draws 1/(1+L) loci
%
%           ConstMag    vector of desired magnitude loci [dB] or 'all'
%                       Default: The same vector like in Matlab plot
%
%           ConstPhase  vector of desired pahse loci [degrees] or 'all'
%                       Default: The same vector like in Matlab plot
%
%   Parameter:
% 
%           color       RGB color triplate; Default: [0.94 0.94 0.94] 
%
%           linestyle   string of desired line style;  Default: ':'
%
%           linewidth   desired line width; Default: 0.8
%
% This file is part of the OpenQsyn toolox, distribted under GNU LGPL3 license

p = inputParser;
p.addOptional('type','T',@(x) validatestring(x,{'T','S'},'qngrid','Type'));
p.addOptional('g_const',[],@(x) validateattributes(x,{'numeric'},{'vector','>',0,'real'},'qngid','ConstMag'));
p.addOptional('p_const',[],@(x) validateattributes(x,{'numeric'},{'vector','real'},'qngid','ConstPhase'));
p.addParameter('color',[0.55 0.55 0.55],@(x) validateattributes(x,{'numeric'},{'size',[1 3]},'qngid','Color'));
p.addParameter('linestyle',':',@(x) validateattributes(x,{'string','char'},{},'qngid','linestyle'));
p.addParameter('linewidth',0.8,@(x) validateattributes(x,{'numeric'},{'scalar','>',0},'qngid','linewidth'));
p.parse(varargin{:});
type = p.Results.type;
g_const = p.Results.g_const;
p_const = p.Results.p_const;
color = p.Results.color;
linestyle = p.Results.linestyle;
linewidth = p.Results.linewidth;

% Figure Handle
figureHandle = gca;
if isempty(figureHandle.Children)
    axObjs.XData(1) = -inf;
    flag = 0;
else
    axObjs = figureHandle.Children;
    flag = 1;
end
if flag
    Ylim(1) = max(figureHandle.YLim(1),-160);
    if Ylim(1) > -20
        Ylim(1) = -20;
    end
    Ylim(2) = min(figureHandle.YLim(2),200);
    if Ylim(2) < 5
        Ylim(2) = 40;
    end
    Xlim = figureHandle.XLim;
    Xlim(1) = ceil(abs(Xlim(1))/45)*45*sign(Xlim(1));
    Xlim(2) = floor(abs(Xlim(2))/45)*45*sign(Xlim(2));
    %if ~rem(Xlim(1)/90,1)
    %    Xlim(1) = figureHandle.XLim(1);
    %elseif -450 < Xlim(1) && Xlim(1) <= 0 % Xlim(1) is between [-450 0]
    %    Xlim(1) = -360; % Xlim(1) is -360 unless if
    %    if (Xlim(2)-Xlim(1)) < 5 % Case of pure integrators
    %        temp = floor(Xlim(1)/360);
    %        Xlim(1) = -360+(temp+1)*180;
    %    end
    %elseif Xlim(1) <= -450 % Xlim(1) is less that -450
    %    temp = floor((Xlim(1)+450)/180);
    %    if mod(temp+ceil((Xlim(1)+450)/180),2)
    %        Xlim(1) = -360+temp*180;
    %    else
    %        Xlim(1) = Xlim(1)-90;
    %    end
    %elseif Xlim(1) > 0
    %    Xlim(1) = 0;
    %end
    n = ceil((Xlim(2)-Xlim(1))/360); % Multiples of 360
    %if n == 1
    %    Xlim(2) = Xlim(1)+360;
    %elseif ~rem(Xlim(1)/90,1)
    %    Xlim(2) = figureHandle.XLim(2);
    %else
    %    Xlim(2) = 180*ceil(Xlim(2)/90);
    %end
    %if Xlim(2) < axObjs(end).XData(1)
    %    Xlim(2) = Xlim(2)+round(Xlim(2)/(axObjs(end).XData(1)+eps))*180;
    %end
    if n == 1 && Xlim(1) < -360
        n = 2;
    end
else
    switch type
        case 'T'
            Ylim = [-100 40];
        case 'S'
            Ylim = [-40 40];
    end
    Xlim = [-360 0];
    n = 1;
    figureHandle.XLim = Xlim;
    figureHandle.YLim = Ylim;
end
% Check if inputs are exists
if isempty(g_const) %strcmp(g_const,'all')
    g_const = [10*(round(Ylim(1)/10)):20:-20 -12 -6 -3 -1 0 0.25 0.5 1 3 6]';
end
if isempty(p_const) %strcmp(p_const,'all')
    p_const = [-180 -150 -120 -90 -50 -30.3 -20.2 -10.1 -5.05 -1.01]';
    p_const = [p_const ; -flip(p_const)];
    p_const(p_const == 180) = [];
end

% Gain Constant Circles around -180 degrees
[Mp,Mg] = NicholsMNcircles(type,g_const,'M');
% Phase Constant Circles around -180 degrees
[Np,Ng] = NicholsMNcircles(type,p_const,'N');
% Critical Point around -180 degrees
CriticalPoint = [-180 0]';

m = 8; % Multiplies the circles by m
Mp = repmat(Mp,1,m);
Mg = repmat(Mg,1,m);
Np = repmat(Np,1,m);
Ng = repmat(Ng,1,m);
CriticalPoint = repmat(CriticalPoint,1,m);
for i = 1:m % create grid from +/- -180-round(m/2)*360
    CriticalPoint(1,i) = -180+(i-round(m/2))*360;
    Mp(:,(i-1)*length(g_const)+1:i*length(g_const)) = Mp(:,(i-1)*length(g_const)+1:i*length(g_const))+180+CriticalPoint(1,i);
    Np(:,(i-1)*length(p_const)+1:i*length(p_const)) = Np(:,(i-1)*length(p_const)+1:i*length(p_const))+180+CriticalPoint(1,i);
end

hold on
plot(Mp,Mg,'color',color,'linestyle',linestyle,'linewidth',linewidth);
plot(Np,Ng,'color',color,'linestyle',linestyle,'linewidth',linewidth);
plot(CriticalPoint(1,:),CriticalPoint(2,:),'+r','linewidth',1);
box on

Ticks = 45*n;
if n > 2
    Ticks = 180;
end

set(figureHandle,'XLim',Xlim)
set(figureHandle,'YLim',Ylim)
set(figureHandle,'XTick',Xlim(1):Ticks:Xlim(2));

switch type
    case 'T'
        Ind = length(Mp);
    case 'S'
        Ind = 1;
end
for i = 1:length(g_const) % Mnichols Text
    for j = 1:m
        if Xlim(1) <= Mp(end,(j-1)*length(g_const)+i) && Mp(end,(j-1)*length(g_const)+i) <= Xlim(2)
            if g_const(i) > 0
                h2 = text(Mp(end,(j-1)*length(g_const)+i),Mg(end,(j-1)*length(g_const)+i),[num2str(g_const(i)),' dB']);
            elseif g_const(i) < 0
                if n < 2
                    h2 = text(Mp(Ind,(j-1)*length(g_const)+i)-35,Mg(Ind,(j-1)*length(g_const)+i)+1,[num2str(g_const(i)),' dB']);
                else
                    h2 = text(Xlim(2)-45,Mg(Ind,(j-1)*length(g_const)+i)+1,[num2str(g_const(i)),' dB']);
                end
            else
                h2 = text(Mp(Ind,(j-1)*length(g_const)+i)+4,Mg(Ind,(j-1)*length(g_const)+i)+1,[num2str(g_const(i)),' dB']);
            end
            set(h2,'fontname','Arial')
            set(h2,'fontsize',9)
        end
    end
end
title('Nichols Chart','fontname','Roman','fontsize',12)
xlabel('Open-Loop Phase (deg)','fontname','Arial','fontsize',10)
ylabel('Open-Loop Gain (dB)','fontname','Arial','fontsize',10)
end

function [Mp,Mg] = NicholsMNcircles(type,pg_const,CircleType)

switch CircleType
    case 'M'
        p_vec = (1.01*pi/180:pi/1000:358.99*pi/180)';
        z = zeros(length(p_vec),length(pg_const));
        OpenLoopCircle = zeros(length(p_vec),length(pg_const));
        Mp = zeros(length(p_vec),length(pg_const));
        for i = 1:length(pg_const)
            z(:,i) = 10.^(pg_const(i)/20)*exp(1i*p_vec);
            switch type
                case 'T'
                    OpenLoopCircle(:,i) = z(:,i)./(1-z(:,i));
                    Mp(:,i) = unwrap(angle(OpenLoopCircle(:,i)))*180/pi-360;
                case 'S'
                    OpenLoopCircle(:,i) = (1-z(:,i))./z(:,i);
                    Mp(:,i) = unwrap(angle(OpenLoopCircle(:,i)))*180/pi;
            end
        end
        Mg = 20*log10(abs(OpenLoopCircle));
    case 'N'
        figureHandle = gca;
        g_vec = (figureHandle.YLim(1)+0.1:0.01:6)';
        z = zeros(length(g_vec),length(pg_const));
        OpenLoopCircle = zeros(length(g_vec),length(pg_const));
        Mp = zeros(length(g_vec),length(pg_const));
        for i = 1:length(pg_const)
            z(:,i) = 10.^(g_vec/20).*exp(1i*pg_const(i)*pi/180);
            switch type
                case 'T'
                    OpenLoopCircle(:,i) = z(:,i)./(1-z(:,i));
                    if imag(z(end,i)) > 0
                        Mp(:,i) = unwrap(angle(OpenLoopCircle(:,i)))*180/pi-360;
                    else
                        Mp(:,i) = unwrap(angle(OpenLoopCircle(:,i)))*180/pi;
                    end
                case 'S'
                    OpenLoopCircle(:,i) = (1-z(:,i))./z(:,i);
                    if imag(z(end,i)) > 0
                        Mp(:,i) = unwrap(angle(OpenLoopCircle(:,i)))*180/pi;
                    else
                        Mp(:,i) = unwrap(angle(OpenLoopCircle(:,i)))*180/pi-360;
                    end
            end
        end
        Mg = 20*log10(abs(OpenLoopCircle));      
end
end