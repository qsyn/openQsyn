function ha = bodeplotter(tpl,w,opt,col,linespec)
%BODEPLOTTER does the main bode ploting
%Subroutine to qtpl.bode and qplant.bodcases
%
%   Inputs
% 
%   TPL         the template points 
%   w           frequencies 
%   opt         how to plot
%   col         colors RGB array
%   linspec     line plotting specs

if nargin<5
    linespec = struct('width',1,'style','-'); % def linespace
end
N = size(tpl,1);

% magnitude
if strcmp(opt,'magphase')
    ha(1)=subplot(2,1,1); 
end

if ~strcmp(opt,'phase')
    for k=1:N
        semilogx(w,imag(tpl(k,:)),'Color',col(k,:),'Tag',num2str(k),...
            'linewidth',linespec.width,'linestyle',linespec.style);
        hold on
    end
    %if strcmp(CASE,'all') || strcmp(CASE,'nom')
    %    [w_nom,t_nom]=gettpl(tplf,'nom');
    %    semilogx(w_nom,imag(t_nom),'--k','linewidth',2);
    %end
    xlim([w(1) w(end)])
    ylabel('Mag [db]')
    
end

% phase
if strcmp(opt,'magphase')
    ha(2)=subplot(2,1,2);
end

if ~strcmp(opt,'mag')
    for k=1:N
        phase=real(tpl(k,:));
        %phase=unwrap(phase*pi/180)*180/pi;
        %if phase(1) > 5
        %    phase=phase-360;
        %end
        semilogx(w,phase,'Color',col(k,:),'Tag',num2str(k),...
             'linewidth',linespec.width,'linestyle',linespec.style); 
         hold on
    end
    %if strcmp(CASE,'all') || strcmp(CASE,'nom')
    %    phase=real(t_nom);
    %    phase=unwrap(phase*pi/180)*180/pi;
    %    if phase(1)>5, phase=phase-360; end
    %    semilogx(w_nom,real(t_nom),'--k','linewidth',2);
    %end
    xlim([w(1) w(end)])
    ylabel('Phase [deg]')
end

xlabel('Frequency [rad/s]')

dcm = datacursormode(gcf);
set(dcm,'updatefcn',@datatipfunc)
end

% custom data cursor tip 
function output_txt = datatipfunc(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
dts = get(event_obj.Target,'Tag');
output_txt = {['Freq: ',num2str(pos(1),4), '[rad/s]'],...
              ['Mag: ',num2str(pos(2),4),' [db]'],...
              ['par set: ',dts]
              };
          
end

