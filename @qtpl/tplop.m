function [ T ] = tplop( A,B,op )
%TPLOP operations between templates 
%
%   Using tplop is not reccomanded and therefore hiiden from the user! 
%   Use tplshift for moving templates across the Nichols chart. 
%
%	[ T ] = tplop( A,B,OP ) performs the operation described by OP between
%	qtpl object A, and an object B.
%
%   Note that the plus opperation is performed in Nichols form (deg+i*db),
%   i.e., for siso transfer functions A,B: plus(A,B) = A(s)*B(s).
%   
%   Exmaple: 
%   [ T ] = tplop( A,B,'+' )   performs an addition between qtpl objects A
%   an dB, returns output as a qtpl object T.


if isa(A,'qtpl')
    tplA = A;
    pnameA = A(1).parNames;
    pnameB ={};
    AisQtpl=1;
else
    % in case A is not a qtpl the second output is. A is copied into B for
    % type checking. In computation loop AisQtpl cariable is used to switch
    % back A and B to tpl1 and tpl2
    tplA = B;
    pnameB = B(1).parNames;
    pnameA = {};
    B = A;
    AisQtpl=0;
end
w = [tplA.frequency]';

N = length(w);

if isnumeric(B) && isscalar(B)
    t2.template = B;
    t2.parameters =[];
    tplB = repmat(t2,N,1);
elseif isa(B,'qfr') || isa(B,'lti')
    t2c = squeeze(freqresp(B,w));
    t2n = c2n(t2c,'unwrap');
    tplB = struct('template',[],'parameters',[]);
    for k=1:N
        tplB(k).template = t2n(k);
        tplB(k).parameters = [];
    end
elseif isa(B,'qtpl')
    if ~all( w==[B.frequency]' )
        error('frequencies must be identical for plus operation')
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
    
    t1 = repmat(tpl1(k).template,length(tpl2(k).template),1);
    t2 = repmat(tpl2(k).template,length(tpl1(k).template),1);
        
    switch op
        case isa(op,'function_handle'),  t = op(t1,t2); 
        case '+', t = t1 + t2;
        case '-', t = t1 - t2;
        case '*', t = t1.*t2;
        case '/', t = t1./t2;
        otherwise, error('undefined operation'); 
    end
    
    if length(tpl2(k).template) > 1
        idx=boundary(real(t),imag(t),0.3);
    else
        idx = 1:length(tpl1(k).template);
    end
    
    n1 = max(1,length(tpl1(k).parameters));     % exclude 0 length
    n2 = max(1,length(tpl2(k).parameters));     % exclude 0 length
    p1 =  repmat(tpl1(k).parameters,n2,1);
    p2 =  repmat(tpl2(k).parameters,n1,1);
    p = [p1 ; p2];
    
    T(k).parameters = p(:,idx);
    T(k).template = t(idx);
    
        
    % if isempty(pnameA) && isempty(pnameB)
    %     T(k).parNames = [];
    % elseif isempty(tpl1(k).parNames)
    %     T(k).parNames = pnameB;
    % elseif isempty(pnameB)
    %     T(k).parNames = tpl1(k).parNames;
    % else
    %     T(k).parNames = { tpl1(k).parNames{:} pnameB{:} };
    % end
    T(k).parNames = [pnameA pnameB];
    
end
        
end