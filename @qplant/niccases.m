function niccases(obj,varargin)
%NICCASES plant frequency response for selected cases on Nichols chart
%   QPLANT.NICCASES replace CASES in Qsyn
%
%   NICCASES(qplant)   plots Nichols chart for all cases given by the
%   plant parameters
%
%   NICCASES(qplant,par,w)   specify frequency and cases set
%
%   NICCASES(qplant,par,w,parameter,value)   specify additional options
%   using parameter/value pairs
%   
%   qplant.NICCASES(...)    alternative usage
%
%
%   Inputs (Optional):
%       par     array with each column a different parameter case;
%               default is the eniter grid
%       w       frequency vector. default is the frequency vector of 
%               nominal case, if available, or logspace(-2,2,50)
%
%   Additional Options:
%       color       RGB array
%       shownom     show nominal in bold dashed line 0 (def) | 1
%
%   See also: qplant/bodcases

p = inputParser;
addOptional(p,'par',[],@(x) validateattributes(x,{'numeric'},{'2d'},1))
addOptional(p,'w',[],@(x) validateattributes(x,{'numeric'},{'positive'},1))
addParameter(p,'color',[],@(x) validateattributes(x,{'numeric'},{'ncols', 3}))
addParameter(p,'shownom',0,@(x) validateattributes(x,{'numeric'},{'scalar','binary'}))
parse(p,varargin{:})

par = p.Results.par;
w = p.Results.w;
col = p.Results.color;
shownom = p.Results.shownom;

[res,w] = cases(obj,par,w);     % compute cases
N = size(res,2);
if isempty(col)                 % pick colors
    col = lines(size(res,2));
elseif size(col,1) < N
    col = repmat(col,ceil(N/size(col,1)),1);
end


linespec = struct('width',1,'style','-');
nicholsplotter(res.',-w,col,linespec);  % plot the Bode for all cases.
                                        % the (-w) is to get the frequency
                                        % data tip dispalyed right.

if  shownom
    hold on
    linespec.width = 3;
    linespec.style = '--';
    nompar = obj.pars.nom;
    nomres = cases(obj,nompar,w);
    nicholsplotter(nomres.',w,[0 0 0],linespec);   % plot the Nichols for nominal case  
end

end