% specs

rsrs('ex1',[],[1.2 0.2],10,1.5,[],logspace(-1,2),2.85,3.1);

add2spc('ex1','odsrs',logspace(-1,2),6)

tic
cbnd('ex1','odsrs');
toc
%%
t1 = P.templates(1);
tpl = t1.template;
gphase=-360:10:0;
gmag=-50:5:50;
specfunc=@fodsrs;
%spec=6;
spec1 = qspc('odsrs',w,6)
spec2 = qspc.rsrs([1.2 0.2],10,1.5,[],logspace(-1,2),2.85,3.1)

des = qdesign(P,[spec1 spec2]) 
%bnd = des.makebnd(tpl,specfunc,6,gphase,gmag)
tic
des.cbnd('odsrs')
toc
%scatter(real(bnd.c{1}),imag(bnd.c{1}))
des.cbnd('rsrs')
