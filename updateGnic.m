function [] = updateGnic(G)
des = evalin('base','des');
Lnom = series(des.nom,G);
plot(app.nicholsPlot,Lnom.phase,Lnom.mag,'k','LineWidth',1.2);
hold(app.nicholsPlot,"on")
for k=1:length(des.tpl)
    str = sprintf(' %g',des.tpl(k).frequency);
    tk = qfr(des.tpl(k).template(1),des.tpl(k).frequency);
    Ltpl = series(tk,G);
    x_text = Ltpl.phase;
    y_text = Ltpl.mag;
    plot(app.nicholsPlot,Ltpl.phase,Ltpl.mag,'marker','square',...
        'markeredgecolor','k','markerfacecolor',t_color(k,:));
    text(app.nicholsPlot,x_text,y_text,str,'clipping','on') % single space added
end
hold(app.nicholsPlot,"off")
end