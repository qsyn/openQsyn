classdef qexpression
    %QFTEXPRESSION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        expression
        pars
    end
    
    methods
        function exp = qexpression(A,B,op)
            if ischar(A) && isa(B,'qpar')
                exp.expression = A;
                exp.pars = B;
            elseif (isa(A,'qpar')) && (isa(B,'qpar'))
                exp.expression = sprintf('%s %s %s',A.name,op,B.name);
                exp.pars = unique(vertcat(A, B));
            elseif  (isa(A,'qexpression')) && (isa(B,'qpar'))
                exp.expression =  sprintf('(%s) %s %s',A.expression,op,B.name);
                exp.pars = unique(vertcat(A.pars, B));
            elseif (isa(A,'qpar')) && (isa(B,'qexpression'))
                exp.expression = sprintf('%s %s (%s)',A.name,op,B.expression);
                exp.pars = unique(vertcat(B.pars, A));
            elseif isnumeric(A) && (isa(B,'qpar'))
                exp.expression = sprintf('%g %s %s',A,op,B.name);
                exp.pars = B;
            elseif isnumeric(B) && (isa(A,'qpar'))
                exp.expression = sprintf('%s %s %g',A.name,op,B);
                exp.pars = A;
            elseif isnumeric(A) && (isa(B,'qexpression'))
                exp.expression = sprintf('%g %s (%s)',A,op,B.expression);
                exp.pars = B.pars;
            elseif isnumeric(B) && (isa(A,'qexpression'))
                exp.expression = sprintf('(%s) %s %g',A.expression,op,B);
                exp.pars = A.pars;
            elseif (isa(A,'qexpression')) && (isa(B,'qexpression'))
                exp.expression = sprintf('(%s) %s (%s)',A.expression,op,B.expression);
                exp.pars = unique(vertcat(A.pars,B.pars));
            end
        end
        function exp = plus(A,B)
            exp = qexpression(A,B,'+');
        end
        function exp = uplus(A)
            exp = qexpression([],A,'+');
        end
        function exp = minus(A,B)
            exp = qexpression(A,B,'-');
        end
        function exp = uminus(A)
            exp = qexpression([],A,'-');
        end
        function obj = mtimes(A,B)
            if isa(A,'lti')
                obj = qftlti(B,A);
            elseif isa(B,'lti')
                obj = qftlti(A,B);
            else
                obj = qexpression(A,B,'*');
            end
        end
        function exp = mrdivide(A,B)
            exp = qexpression(A,B,'/');
        end
        function h = qexp2func(obj)
            args = sprintf('%s, ',obj.pars.name);
            argF = sprintf('@(%s)  ',args(1:end-2));
            %exp = replace(obj.expression,{'*','/','^'},{'.*','./','.^'}); %introduced on 2016b
            exp = strrep(strrep(strrep(obj.expression, '*', '.*'),'/','./'),'^','.^');
            h = str2func([argF exp]);
        end
        function val = nom(obj)
            %try 
            %    se = str2sym(obj.expression);
            %    val = double(subs(se,{obj.pars.name},{obj.pars.nominal}));
            %catch
            %   disp('no symbolic toolbox avilable --> using the evil STR2FUNC')
            fh =  qexp2func(obj);
            val = fh(obj.pars.nominal);
            %end
        end
        function val = expcase(obj,par)
            %EXPCASE return the QEXPRESSION value for a given parameter
            %case
            fh =  qexp2func(obj);
            PAR = num2cell(par);
            val = fh(PAR{:});
            %end
        end
        function obj = horzcat(varargin)
            obj = qpoly(varargin{:});
        end
    end
end

