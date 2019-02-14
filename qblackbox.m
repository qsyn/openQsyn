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
    end
    
       methods(Static)
        function T = funcval(f,c,s)
            %FUNCVAL return the value of the plant, represented by a
            %function handle F created by QPLANT2FUNC, for given parameter
            %case C and for S. The output is in unwrapped Nichols format.
            ck = c;
            ck{end} = s;
            T = zeros(length(s),length(c{1}));
            pargrid = cell2mat(c(1:end-1));
            for k=1:size(pargrid,2) % for each parameter set...
                ck(1:size(pargrid,1)) = mat2cell(pargrid(:,k),ones(1,size(pargrid,1)));
                nyq  = f(ck{:});
                T(:,k) = c2n(nyq,'unwrap');
            end
        end
    end
 
    
end

