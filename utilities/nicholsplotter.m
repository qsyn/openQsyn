function ha = nicholsplotter(tpl,w,col,linespec)
%NICHOLSPLOTTER does the main Nichols ploting
%Subroutine to qtpl.niccases and
%
%   Inputs
% 
%   TPL     the template points 
%   w       frequencies 
%   col     colors RGB array


N = size(tpl,1);
for k=1:N
    plot(real(tpl(k,:)),imag(tpl(k,:)),'Zdata',w,...
        'Color',col(k,:),'Tag',num2str(k),...
        'linewidth',linespec.width,'linestyle',linespec.style);
    hold on
end
hold off   
ha = gca;
xlabel('Phase [deg]')
ylabel('Magnitude [dB]')

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
output_txt = {['Phase: ',num2str(pos(1),4), '[deg]'],...
              ['Mag: ',num2str(pos(2),4),' [dB]'],...
              ['frequency: ',num2str(abs(pos(3)),4),' [rad/s]'],...
              ['par set: ',dts]
              };
          
end

