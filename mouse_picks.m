function [h] = mouse_picks(G,Controller)
fprintf('%d',Controller.gain(1));
h = figure;
[p,z] = pzmap(G);
plot(real(p),imag(p),'xr','buttondownfcn',{@Mouse_Callback_p,'down'});
hold on
plot(real(z),imag(z),'ob','buttondownfcn',{@Mouse_Callback_z,'down'});
hold on
plot([-5 5],[0 0],'k','linewidth',1)
hold on
plot([0 0],[-5 5],'k','linewidth',1)
grid minor
hold on
text(-2.5,2.5,'Close figure after changing poles / zeroes')

% Callback function for poles
function Mouse_Callback_p(hObj,~,action)

persistent curobj xdata ydata ind 
pos = get(gca,'CurrentPoint');

switch action
    case 'down'
        curobj = hObj;
        xdata = get(hObj,'xdata');
        ydata = get(hObj,'ydata');
        [~,ind] = min(sum((xdata-pos(1)).^2+(ydata-pos(3)).^2,1));
        set(gcf,...
            'WindowButtonMotionFcn',  {@Mouse_Callback_p,'move'},...
            'WindowButtonUpFcn',      {@Mouse_Callback_p,'up'});
        
    case 'move'
        if ydata(ind) ~= 0
            if ~mod(ind,2) == 1
                ind = ind - 1;
            end
            ydata(ind+1) = - ydata(ind);
            xdata(ind+1) = xdata(ind);
            ydata(ind) = pos(3);
            set(curobj,'ydata',ydata)
            xdata(ind) = pos(1);
            set(curobj,'xdata',xdata)
        else
            xdata(ind) = pos(1);
            set(curobj,'xdata',xdata)
        end
    case 'up'
        set(gcf,...
            'WindowButtonMotionFcn',  '',...
            'WindowButtonUpFcn',      '');
end
if ~isempty(ind)
    assignin("base","newPos",{[pos(1) pos(3)] ind 'pole'})
end
end

% Callback function for zeros
function Mouse_Callback_z(hObj,~,action)

persistent curobj xdata ydata ind 
pos = get(gca,'CurrentPoint');

switch action
    case 'down'
        curobj = hObj;
        xdata = get(hObj,'xdata');
        ydata = get(hObj,'ydata');
        [~,ind] = min(sum((xdata-pos(1)).^2+(ydata-pos(3)).^2,1));
        set(gcf,...
            'WindowButtonMotionFcn',  {@Mouse_Callback_z,'move'},...
            'WindowButtonUpFcn',      {@Mouse_Callback_z,'up'});
        
    case 'move'
        if ydata(ind) ~= 0
            if ~mod(ind,2) == 1
                ind = ind - 1;
            end
            ydata(ind+1) = - ydata(ind);
            xdata(ind+1) = xdata(ind);
            ydata(ind) = pos(3);
            set(curobj,'ydata',ydata)
            xdata(ind) = pos(1);
            set(curobj,'xdata',xdata)
        else
            xdata(ind) = pos(1);
            set(curobj,'xdata',xdata)
        end
    case 'up'
        set(gcf,...
            'WindowButtonMotionFcn',  '',...
            'WindowButtonUpFcn',      '');
end
if ~isempty(ind)
    assignin("base","newPos",{[pos(1) pos(3)] ind 'zero'})
end
end

end
