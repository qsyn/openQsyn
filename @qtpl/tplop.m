function [ T ] = tplop( A,B,op )
%TPLOP operations between templates 
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

tpl1 = A;
w = [tpl1.frequency]';
N = length(w);

if isnumeric(B) && isscalar(B)
    t2.template = B;
    t2.parameters =[];
    tpl2 = repmat(t2,N,1);
elseif isa(B,'qfr') || isa(B,'lti')
    t2c = squeeze(freqresp(B,w));
    t2n = c2n(t2c,'unwrap');
    tpl2 = struct('template',[],'parameters',[]);
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
    
    p1 =  repmat(tpl1(k).parameters,length(tpl2(k).parameters),1);
    p2 =  repmat(tpl2(k).parameters,length(tpl1(k).parameters),1);
    p = [p1 p2];
    
    T(k).parameters = p(:,idx);
    T(k).template = t(idx);
    
end
        
end