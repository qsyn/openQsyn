function [ T ] = cpop( A,B,opr,varargin )
%CPOP operations between templates in the complex plain
%
%	[ T ] = CPOP( A,B,OPR ) performs the operation described by OPR between
%	qtpl object A, and an object B.
%
%   Note that the plus opperation is performed in complex form (re+im*i),
%   thus templates are transformed to complex form and back
%
%   Inputs 
%   A       qtpl array
%   B       one of the following: 
%               qtpl array with same frequenciers as A
%               numeric array of frequency response points in complex form
%               Matlab control toolboxz LTI siso object
%               qfr object               
%   opr     operation: '+'|'-'|'*'|'/'|'sens'|'comp'
%
%   Output
%   T       qplt array
%   
%   Exmaple: 
%   [ T ] = cpop( A,B,'+')   performs an addition between qtpl objects A
%   and B, returns output as the qtpl object T.

%-- Parse inputs --
p = inputParser;
addRequired(p,'A',@(x)isa(x,'qtpl'));
addRequired(p,'B');
addRequired(p,'opr',@(x)validateattributes(x,{'char'},{'nonempty'}));
addParameter(p,'prunning',0,@(x)validateattributes(x,...
    {'numeric'},{'>=',0,'<=',2,'integer','scalar'},'cpop','option')); % set prunning method
parse(p,A,B,opr,varargin{:});

tpl1 = p.Results.A;
B = p.Results.B;
opr = p.Results.opr;
prunning = p.Results.prunning;
%-----------------

w = [tpl1.frequency]';
N = length(w);

if isnumeric(B)
    t2.template = c2n(B);
    t2.parameters =[];
    tpl2 = repmat(t2,N,1);
elseif isa(B,'qfr') || isa(B,'lti')
    t2c = squeeze(freqresp(B,w));
    t2n = c2n(t2c,'unwrap');
    tpl2 = struct('template',t2n,'parameters',[]);
    for k=1:N
        tpl2(k).template = t2n(k);
        tpl2(k).parameters = [];
    end
elseif isa(B,'qtpl')
    if ~all( w==[B.frequency]' )
        error('frequencies must be identical for plus operation')
    end
    tpl2 = B;
else
    error(['second argument must be either a numeric scalar, ',... 
        'QTPL object, QFR object, or LTI object']);
end

T =qtpl(N);
for k = 1:N
    
    T(k).frequency = w(k);
    
    t1c = n2c(tpl1(k).template);
    t2c = n2c(tpl2(k).template);
    n1 = length(t1c);
    n2 = length(t2c);
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
    
    if length(tpl2(k).template) > 1 && prunning ==1
        [~,idx] = prune(t,[2 2]);
    elseif length(tpl2(k).template) > 1 && prunning == 2
        idx=boundary(real(t),imag(t),0.3);
    else
        idx = 1:length(t1);
    end
    
    p1 =  repmat(tpl1(k).parameters,1,n2);
    p2 =  repmat(tpl2(k).parameters,1,n1);
    p = [p1 p2];
    
    T(k).parameters = p(:,idx);
    T(k).template = t(idx);
    
end
        
end