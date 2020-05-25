function newPos = mouse_picks(G,Controller)
figure(1)
[p,z] = pzmap(G);

plot(real(p),imag(p),'xr','buttondownfcn',{@Mouse_Callback,'down'});
hold on
plot(real(z),imag(z),'ob','buttondownfcn',{@Mouse_Callback,'down'});
hold on
plot([-10 10],[0 0],'k','linewidth',1)
hold on
plot([0 0],[-10 10],'k','linewidth',1)
grid minor
hold on


% Callback function for each point
function Mouse_Callback(hObj,~,action)

persistent curobj xdata ydata ind
pos = get(gca,'CurrentPoint');

switch action
    case 'down'
        curobj = hObj;
        xdata = get(hObj,'xdata');
        ydata = get(hObj,'ydata');
        [~,ind] = min(sum((xdata-pos(1)).^2+(ydata-pos(3)).^2,1));
        set(gcf,...
            'WindowButtonMotionFcn',  {@Mouse_Callback,'move'},...
            'WindowButtonUpFcn',      {@Mouse_Callback,'up'});
        
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
newPos = pos;
end
end
