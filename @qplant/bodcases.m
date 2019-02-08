function bodcases(obj,varargin)
%BODCASES plant frequency response for selected cases on Bode plot
%   Qplant.BODCASES replaces CASES in Qsyn
%
%   BODCASES(plant)   plots Bode for all cases given by the
%   plant parameters
%
%   BODCASES(plant,par,w)   specify frequency and cases set
%
%   BODCASES(plant,par,w,parameter,value)   specify additional options
%   using parameter/value pairs
%
%   qplant.BODCASES(...)    alternative usage
%
%
%   Inputs (Optional):
%       par     array with each column a different parameter case;
%               default is the eniter grid
%       w       frequency vector. default is the frequency vector of 
%               nominal case, if available, or logspace(-2,2,50)
%
%   Additional Options:
%       color       specifiy color: RGB array
%       shownom     show nominal in bold dashed line: 0 (def) | 1
%       showmag     show magnitude plot: 0 | 1 (def)
%       shophase    show phase plot: 0 | 1 (def)
%
%   See also: qplant/niccases

% parse inputs
p = inputParser;
addOptional(p,'par',[],@(x) validateattributes(x,{'numeric'},{'2d'},1))
addOptional(p,'w',[],@(x) validateattributes(x,{'numeric'},{'positive'},1))
addParameter(p,'color',[],@(x) validateattributes(x,{'numeric'},{'ncols', 3}))
addParameter(p,'showmag',1,@(x) validateattributes(x,{'numeric'},{'scalar','binary'}))
addParameter(p,'showphase',1,@(x) validateattributes(x,{'numeric'},{'scalar','binary'}))
addParameter(p,'shownom',0,@(x) validateattributes(x,{'numeric'},{'scalar','binary'}))
parse(p,varargin{:})

par = p.Results.par;
w = p.Results.w;
col = p.Results.color;
shownom = p.Results.shownom;

if p.Results.showmag && p.Results.showphase
    opt = 'magphase';
elseif p.Results.showmag
    opt = 'mag';
else
    opt = 'phase';
end

[res,w] = cases(obj,par,w);     % compute cases
N = size(res,2);
if isempty(col)                 % pick colors
    col = lines(size(res,2));
elseif size(col,1) < N
    col = repmat(col,ceil(N/size(col,1)),1);
end

linespec = struct('width',1,'style','-');
bodeplotter(res.',w,opt,col,linespec);   % plot the Bode for all cases

if  shownom
    linespec.width = 3;
    linespec.style = '--';
    nompar = obj.pars.nom;
    nomres = cases(obj,nompar,w);
    bodeplotter(nomres.',w,opt,[0 0 0],linespec); % plot the Bode for nominal case
end

end
