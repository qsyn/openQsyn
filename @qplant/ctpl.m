function obj = ctpl(obj,method,w,varargin)
%CTPL computes the templates for given qplant object
%
%   P = ctpl(P,method,w)    compute templates for plant P vai specified
%   method at given frequecies.
%
%	P = ctpl(P,method,w,parameter,value)    specify additional options
%	using parameter/value pairs
%
%   Inputs (required):
%
%       obj     qplant object
%
%       method  template computation method. one of the following:
%               'grid'      uniform grid over parameter
%               'rngrid'    random grid over the parameters
%               'random'    random sample points
%               'recgrid'   recurcive grid (TO DO)
%               'aedgrid'   recurcive edge grid
%               'cases'     explicitly given parameter cases
%
%       w       frequency vector in rad/s
%
%   Additional Options:
%
%       plotOn      plots during the computation. 0 (def) | 1.
%       accuracy    set accurcy for computation as [deg_accuracy , dB_accuracy]
%                   for recgrid and aedgrid methods. def = [5 3].
%       union       unite new tpls with existing ones. 0 (def) | 1
%       parameters  parameter set for cases method
%       size        size of grid for methods 'grid', 'rngrid' and 'random' 
%   
%   See also: qplant/cnom

if nargin<4, options=[]; end

% parse inputs
p = inputParser;
addRequired(p,'method',@(x) validateattributes(x,{'char'},{'nonempty'}));
addRequired(p,'w',@(x) validateattributes(x,{'numeric'},{'vector','positive','real'}));
addParameter(p,'plotOn',0,@(x) validateattributes(x,{'numeric'},{'scalar','binary'}));
addParameter(p,'union',0,@(x) validateattributes(x,{'numeric'},{'scalar','binary'}));
addParameter(p,'accuracy',[5 3],@(x) validateattributes(x,{'numeric'},{'size',[1 2],'positive','real'}));
addParameter(p,'parameters',[],@(x) validateattributes(x,{'numeric'},{'2d'}));
addParameter(p,'size',[],@(x) validateattributes(x,{'numeric'},{'2d'}));
parse(p,method,w,varargin{:})

method = p.Results.method;
w = p.Results.w;
options.plot_on = p.Results.plotOn;
options.Tacc = p.Results.accuracy;
options.pars = p.Results.parameters;
N = p.Results.size;

% compute based on given method:
switch method
    case 'grid', tpl=obj.cgrid(w,0,N);
    case 'rndgrid', tpl=obj.cgrid(w,1,N);
    case 'random', tpl=obj.cgrid(w,2,N);
    case 'recgrid', tpl=obj.recgrid(w,options);
    case 'recedge', tpl=obj.recedge(w,options);
    case 'cases', tpl=obj.cases2tpl(options,w);
    otherwise, error('unrecognized method!')
end

% add parameter names and the nominal point (at beginning of each template)
pname = {obj.pars.name};
pnom = [obj.pars.nominal]';
tnom=cases(obj,pnom,w);
for k=1:length(w)
    tpl(k).parNames = pname;                    % add names
    tpl(k)=add2tpl(tpl(k),tnom(k),pnom,'x');    % add nominal point
    tpl(k).unwrap;                              % unwrap according to nominal point
end

% adds uncetin poles/zeros
if any(strcmp({obj.pars.name},'uncint_par'))
    %tplop('A+B',t_,c2n((j*w_tpl(i)).^n_dif(1:(length(n_dif)-1)),'unwrap'));
    idx = strcmp({obj.pars.name},'uncint_par');
    t_uncint = c2n( (1j*w).^obj.pars(idx).discrete );
    tpl = cpop(tpl,t_uncint,'+');
end

% adds unstructured uncertainty
if ~isempty(obj.unstruct)
    tpl = addunstruct(obj,tpl);
end

% if union is on
if p.Results.union 
    tpl = union(tpl,obj.templates);
end

if isempty(obj.templates)
    obj.templates = tpl;
else
    % replacing existing freqeucies and inserting new ones in
    % the right places (sorted by frequency)
    w0 = [obj.templates.frequency];
    w1 = [tpl.frequency];
    w = unique([w0 w1]);
    inew = ismember(w,w1);
    [~,i0] = ismember(w(~inew),w0);
    TPL = qtpl(length(w));
    TPL(inew) = tpl;
    TPL(~inew) = obj.templates(i0);
    obj.templates = TPL;
end


end