function [obj] = add2tpl(obj,tpl,par,opt)
%ADD2TPL adds points to  aqtpl object
% 
%Usage: 
%
%   tplObj = ADD2TPL(oldTplObj,tplPoints,Pars,Option) 
%
%Inputs:
%
%   oldTplObj   old qtpl object we want to add the points to
%   tplPoints   vector template points we wish to add 
%   Pars        matrix of parameters coresponding to the added points  
%   Options     'a' (append to end) | 'x' (insert at start)
%

if (size(tpl,1))~=size(par,2)
    error('PAR must have same number of columns as the raws in TPL')
end


if strcmp(opt,'a')
    obj.template = [obj.template ; tpl];
    obj.parameters = [obj.parameters par];
elseif strcmp(opt,'x')
    obj.template = [tpl ; obj.template];
    obj.parameters = [par obj.parameters];
else
    error('options: ''a'' (append to end) | ''x'' (insert at start)')
end


end

