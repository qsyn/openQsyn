function tpl = cgrid(obj,w,rnd,n)
%CGRID computes tpl for a qsys element using simple grid
%Facilitates grid, random grid, and random sampling
%
%For random and random grid the parameter set is random, yet identical in
%each frequency.
%
% Usage:    
%   tpl = cgrid(T,w,0,n)    computes templates for qsys element T at
%   frequencies specified by w using a uniform grid of n points 
%
%   tpl = cgrid(T,w,1,n)  computes using a random grid of parameters
%
%   tpl = cgrid(T,w,2,n)  computes using randomly sampled parametes 
%  

if nargin<3, rnd=0; end
idx = ~strcmp({obj.pars.name},'uncint_par');
switch rnd
    case 0
        methodName='grid';
        pgrid = grid(obj.pars(idx),n,rnd);
    case 1
        methodName='random grid';
        pgrid = grid(obj.pars(idx),n,rnd);
    case 2
        methodName='random sampling';
        pgrid = sample(obj.pars(idx),n); % correct usage: options.cases(=100)
    otherwise, error('unavilable rnd option')
end

fprintf('Calculating templates using the %s method \n',methodName)

%col = distinguishable_colors(length(w)); % to remove
%col = lines(length(w));
tpl = qtpl(length(w)); % pre-allocating
for k = 1:length(w)
    fprintf('--> for w=%g [rad/s] \n',w(k));
    T = cases(obj,pargrid,w(k));
    tpl(k)=qtpl(w(k),T,pgrid);
    %scatter(real(T),imag(T),10,col(k,:)); hold on
end

end