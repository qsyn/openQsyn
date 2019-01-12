col = distinguishable_colors(length(w));

h = des.showbnd('odsrs',[],[5 10 20 50],col(5:end,:));
des.showbnd('rsrs',h,[0.2 0.5 1 2],col(1:4,:));

s = zpk(0,[],1);
set(s,'DisplayFormat','Frequency')
G = 2.5*(1+s/6)*(1+2*0.6*s/4+s^2/16)/s/(1+s)/(1+s/3.2)/(1+s/26)

des.loopnic(G)

%% prefilter design
T = des.tpl.comp;

spec2.show('freq')
des.clmag(G,1)
ylim([-55 10])

%%
F = 1/(1+2*0.83*s/3.4+s^2/3.4^2)
spec2.show('freq')
des.clmag(G,F)
ylim([-55 10])