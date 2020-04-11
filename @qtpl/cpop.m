function [ T ] = cpop( A,B,opr,varargin )
%CPOP operations between templates in the complex plain
%
%	T = CPOP( A,B,OPR ) performs the operation described by OPR between
%	A and B, which at least one of them is a qtpl object
%
%   Note that the opperation is performed in complex form (re+im*i),
%   thus templates are transformed to complex form and back
%
%   uncertainties are treated as non-dependent! 
%
%   Inputs 
%   A           qtpl array
%   B           one of the following: 
%               * qtpl array with same frequenciers as A
%               * numeric array of frequency response points in complex form
%               * Matlab control toolboxz LTI siso object 
%               * qfr object               
%   opr         operation: '+'|'-'|'*'|'/'|'sens'|'comp'
%   pruning    optional pruning option: 0 (def) | 1  | 2 
%                   option 1 -- use Qsyn prune function
%                   option 2 -- use Matlab boundary function (from 2017b)
%               
%
%   Output
%   T       qplt array
%   
%   Exmaple: 
%   T = CPOP( A,B,'+')   performs an addition between qtpl objects A
%   and B, returns output as the qtpl object T.
%
%   See also: qtpl/plus qtpl/minus qtpl/times qtpl/rdivide qtpl/sens
%   qtpl/comp


%-- Parse inputs --
p = inputParser;
addRequired(p,'A',@(x)isa(x,'qtpl') || isnumeric(x));
addRequired(p,'B');
addRequired(p,'opr',@(x)validateattributes(x,{'char'},{'nonempty'}));
addParameter(p,'pruning',0,@(x)validateattributes(x,...
    {'numeric'},{'>=',0,'<=',2,'integer','scalar'},'cpop','option')); % set prunning method
parse(p,A,B,opr,varargin{:});

if isa(p.Results.A,'qtpl')
    tplA = p.Results.A;
    pnameA = p.Results.A(1).parNames;
    pnameB ={};
    AisQtpl=1;
else
    % in case A is no a qtpl the second output is. A is copied into B for
    % type checking. In computation loop AisQtpl cariable is used to switch
    % back A and B to tpl1 and tpl2
    tplA = p.Results.B;
    pnameB = p.Results.B(1).parNames;
    pnameA = {};
    B = p.Results.A;
    AisQtpl=0;
end

%tpl1 = p.Results.A;
%B = p.Results.B;
opr = p.Results.opr;
pruning = p.Results.pruning;
%-----------------

w = [tplA.frequency]';
N = length(w);

if isnumeric(B)
    t2.template = c2n(B);
    t2.parameters =[];
    tplB = repmat(t2,N,1);
elseif isa(B,'qfr') || isa(B,'lti')
    t2c = squeeze(freqresp(B,w));
    t2n = c2n(t2c,'unwrap');
    tplB = struct('template',t2n,'parameters',[]);
    for k=1:N
        tplB(k).template = t2n(k);
        tplB(k).parameters = [];
    end
elseif isa(B,'qtpl')
    if ~all( w==[B.frequency]' )
        error('frequencies must be identical for complex plain operations')
    end
    tplB = B;
    pnameB = tplB(1).parNames;
else
    error(['second argument must be either a numeric scalar, ',... 
        'QTPL object, QFR object, or LTI object']);
end

if AisQtpl
    tpl1=tplA;
    tpl2=tplB;
else
    tpl1=tplB;
    tpl2=tplA;
end

if ~iscell(pnameA), pnameA={ pnameA }; end
if ~iscell(pnameB), pnameB={ pnameB }; end

T =qtpl(N);
for k = 1:N
    
    T(k).frequency = w(k);
    t1c = n2c(tpl1(k).template);
    t2c = n2c(tpl2(k).template); 
    
    n1 = length(t1c);
    n2 = length(t2c);
    
    p1 =  repmat(tpl1(k).parameters,1,n2);
    p2 =  repmat(tpl2(k).parameters,1,n1);
    p = [p1 ; p2];
   
    t1 = repmat(t1c,n2,1);
    t2 = repmat(t2c,n1,1);
    
    switch opr
        case isa(opr,'function_handle'),  t = opr(t1,t2); 
        case '+', tc = t1 + t2;
        case '-', tc = t1 - t2;
        case '*', tc = t1.*t2;
        case '/', tc = t1./t2;
        case 'sens', tc = 1./(1 + t1.*t2);           
        case 'comp', tc = t1.*t2./(1 + t1.*t2);  % complementry sensitivity
        otherwise, error('undefined operation'); 
    end
    
    t = c2n(tc,'unwrap');
    
    if  pruning == 1 && length(t) > 1  
        [~,idx] = prune(t,[2 2]);
        idx = [1 idx]; % nominal stays 
    elseif pruning == 2 && length(t) > 1  
        idx=boundary(real(t),imag(t),0.3);
        idx = [1; idx]; % nominal stays
    else
        idx = 1:length(t1);
    end
    
    T(k).parameters = p(:,idx);
    T(k).template = t(idx);
    
    T(k).parNames = [pnameA pnameB];
    
end
        
end