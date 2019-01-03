function [obj] = add2tpl(obj,tpl,par,opt)
%ADD2TPL Summary of this function goes here
%   Detailed explanation goes here

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
    error('options: ''a'' (append to end) | ''i'' (insert at start)')
end


end

