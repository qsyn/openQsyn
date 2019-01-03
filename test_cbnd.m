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
spec = qspc('odsrs',w,6)

des = qdesign(P,spec) 
%bnd = des.makebnd(tpl,specfunc,6,gphase,gmag)
tic
des.cbnd('odsrs')
toc
%scatter(real(bnd.c{1}),imag(bnd.c{1}))