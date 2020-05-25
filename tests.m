clear 
clc
G = tf([1 2 4],[ 1 -2 3 9 1 1 1 ]);
[p,z] = pzmap(G);

I = figure()
plot([-5 5],[0 0],'k','linewidth',1)
hold on
plot([0 0],[-5 5],'k','linewidth',1)
grid minor
j = 1;
for i = 1 : length(p)
    h(j) = images.roi.Point(gca,'Position',[real(p(i)) imag(p(i))],'Color','r');
    j = j+1
end
j = 1;
for i = 1 : length(z)
    h(j) = images.roi.Point(gca,'Position',[real(z(i)) imag(z(i))],'Color','b');
    j = j+1
end
