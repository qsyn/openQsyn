classdef qblackbox < qplant
%QBLACKBOX is a subclass of the qplant class 
    
    properties
        blackBox        % function handle for "black box" plants
    end
    
    methods
        function obj = qblackbox(func,pars)
            %QBLACKBOX Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@qplant(func,pars);
        end
        
        function h = qplant2func(obj)
            %QPLANT2FUNC return an handle to a function object f@(p1,p2,...pn,s)
            %with p1,...,pn corresponding to uncertain parameters.
            %
            %Used by: FUNCVAL, CASES, CTPL, CNOM and their subroutines          
            h = obj.blackBox;
        end
        function obj = cnom(obj,w)
            %CPNOM computes the nominal transfer function
            %
            %   CNOM(P,W)   computes the nominal for qplant object P for the frequency
            %   vector W; results is stored under the property 'nominal'.
            %
            %   See also: qplant/ctpl
            if nargin<2, w = logspace(-1,2,50); end
            f = qplant2func(obj);
            pnom = [obj.pars.nominal];
            N = length(w);
            C=cell(length(pnom),1);
            for k=1:length(pnom)
                C{k} = repmat(pnom(k),N,1);                
            end
            C{end+1}=1j*w;
            nyq = f(C{:});
            obj.nominal=qfr(c2n(nyq,'unwarp'),w) ;
        end
    end
    
       methods(Static)
        function T = funcval(f,c,s)
            %FUNCVAL return the value of the plant, represented by a
            %function handle F created by QPLANT2FUNC, for given parameter
            %case C and for S. The output is in unwrapped Nichols format.
            ck = c;
            ck{end} = reshape(s,1,[]);
            T = zeros(length(s),length(c{1}));
            pargrid = cell2mat(c(1:end-1));
            for k=1:size(pargrid,2) % for each parameter set...
                ck(1:size(pargrid,1)) = mat2cell(repmat(pargrid(:,k),1,length(s)),ones(1,size(pargrid,1)));
                nyq  = f(ck{:});
                T(:,k) = c2n(nyq,'unwrap');
            end
        end
    end
 
    
end

