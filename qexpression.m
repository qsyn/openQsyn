classdef qexpression
    %QFTEXPRESSION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        expression
        pars
    end
    
    methods
        function exp = qexpression(A,B,op)
            if (isa(A,'qftpar')) && (isa(B,'qftpar'))
                exp.expression = sprintf('%s %s %s',A.name,op,B.name);
                exp.pars = [A ; B];
            elseif  (isa(A,'qftexpression')) && (isa(B,'qftpar'))
                exp.expression =  sprintf('(%s) %s %s',A.expression,op,B.name);
                exp.pars = unique([A.pars ; B]);
            elseif (isa(A,'qftpar')) && (isa(B,'qftexpression'))
                exp.expression = sprintf('%s %s (%s)',A.name,op,B.expression);
                exp.pars = unique([B.pars ; A]);
            elseif isnumeric(A) && (isa(B,'qftpar'))
                exp.expression = sprintf('%g %s %s',A,op,B.name);
                exp.pars = B;
            elseif isnumeric(B) && (isa(A,'qftpar'))
                exp.expression = sprintf('%s %s %g',A.name,op,B);
                exp.pars = A;
            elseif isnumeric(A) && (isa(B,'qftexpression'))
                exp.expression = sprintf('%g %s (%s)',A,op,B.expression);
                exp.pars = B.pars;
            elseif isnumeric(B) && (isa(A,'qftexpression'))
                exp.expression = sprintf('(%s) %s %g',A.expression,op,B);
                exp.pars = A.pars;
            elseif (isa(A,'qftexpression')) && (isa(B,'qftexpression'))
                exp.expression = sprintf('(%s) %s (%s)',A.expression,op,B.expression);
                exp.pars = unique([A.pars ; B.pars]);
            end
        end
        function exp = plus(A,B)
            exp = qftexpression(A,B,'+');
        end
        function exp = minus(A,B)
            exp = qftexpression(A,B,'-');
        end
        function obj = mtimes(A,B)
            if isa(A,'lti')
                obj = qftlti(B,A);
            elseif isa(B,'lti')
                obj = qftlti(A,B);
            else
                obj = qftexpression(A,B,'*');
            end
        end
        function exp = mrdivide(A,B)
            exp = qftexpression(A,B,'/');
        end
    end
end

