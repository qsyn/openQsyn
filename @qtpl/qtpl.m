classdef qtpl
    %QFTTPL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frequency   (1,1)   {mustBeNumeric,mustBeReal,mustBeNonnegative}
        template    (:,1)   {mustBeNumeric}
        parameters  (:,:)   {mustBeNumeric}
        parNames
    end
    
    methods
        function obj = qtpl(frequency,template,varargin)
            %QTPL construct a qtpl object
            %
            % Usage
            %
            % obj = QTPL(frequency,template)   constructs a qtpl object for
            % given frequency with given points
            %
            % obj = QTPL(frequency,template,parameter)   also specifies parameter
            % corresponding to the template points 
            %
            % obj = QTPL(frequency,template,parameter,parNames)   also specifies  
            % parameter names

            if nargin==0, return; end
            if nargin==1 % build an empty array of length given by input 1
                obj(frequency,1) = obj;
                return
            end
            
            p = inputParser;
            addRequired(p,'frequency',@(x)validateattributes(x,{'numeric'},{'scalar','positive','real'}));
            addRequired(p,'template',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
            addOptional(p,'parameters',[],@(x)validateattributes(x,{'numeric'},{'2d','real'}));
            addOptional(p,'parNames',{},@(x)validateattributes(x,{'cell'},{'vector'}));
            parse(p,frequency,template,varargin{:});
            
            obj.frequency = p.Results.frequency;
            obj.template = p.Results.template;
            obj.parameters = p.Results.parameters; 
            obj.parNames = p.Results.parNames; 
            
            % complete missing parameters
            if isempty(obj.parameters)
                obj.parameters = 1:length(obj.template);
                obj.parNames = 'index (def)';
            end
            
        end 
        function [ T ] = shift( A,b )
           %SHIFT shifts a qtpl object accross the Nichols chart
           %
           %    C = shift(A,b)  shift the qtpl object A by the scalar b,
           %    given in Nichols form deg+j*db. 
           %    
           T = tplop(A,b,'+');
        end
        function [ T ] = plus( A,B )
            %PLUS adds two qtpl arrays
            %
            %   T = A+B  performs an addition between qtpl objects A and
            %   B, returns output as a qtpl object C.
            %
            %   C = plus( A,B )   ALTERNATIVE EXECUTION
            %
            %   Note that the times opperation is performed in complex domain
            %  
            %   See also: qtpl/minus qtpl/times qtpl/rdivide qtpl/cpop
            T = cpop(A,B,'+');
        end
        function [ T ] = minus( A,B )
            %MINUS substruct two qtpl arrays
            %
            %   C = A-B   performs a substraction between qtpl objects A and B, 
            %   returns output as a qtpl object C.
            %
            %   C = minus(A,B)  alternative execution
            % 
            %   Note that the times opperation is performed in complex domain
            %  
            %   See also: qtpl/plus qtpl/times qtpl/rdivide qtpl/cpop
            T = cpop(A,B,'-');
        end
        function [ T ] = uminus( A )
            %UMINUS unary minus for a qtpl array
            %
            %   C = -A   negates the elements of qtpl object A and stores
            %   the reuislts in C
            %
            %   C = uminus( A )   alternative execution
            %
            %   Note that the times opperation is performed in complex domain
            %  
            %   See also: qtpl/plus qtpl/times qtpl/rdivide qtpl/cpop
            T = cpop(0,A,'-');
        end
        function [ T ] = times( A,B )
            %TIMES elemnt-wise multipications of two qtpl arrays (.*)
            %
            %   C=A.*B   computes the product of qtpl objects A and B, returns 
            %   output as a qtpl object C.
            %
            %   C = times(A,B)  alternative execution
            %
            %   Note that the times opperation is performed in complex domain
            %  
            %   See also: qtpl/plus qtpl/minuss qtpl/rdivide qtpl/cpop
            T = cpop(A,B,'*');
        end
        function [ T ] = rdivide( A,B )
            %RDIVIDE elemnt-wise division of two qtpl arrays (./)
            %
            %   C = A./B   divides each element of qtpl object A by the 
            %   corresponding element in qtpl C 
            %
            %   C = rdivide( A,B )  alternative execution
            %
            %   Note that the times opperation is performed in complex domain
            %  
            %   See also: qtpl/plus qtpl/minus qtpl/times qtpl/cpop
            T = cpop(A,B,'/');
        end
        function [ T ] = sens(A,B)
            %SENS compute template of the sensitivinty trnasfer function 
            if nargin<2, B = 1; end
            T = A.cpop(B,'sens');
        end
        function [ T ] = comp(A,B)
            %COMP compute template of the comp. sensitivinty trnasfer function 
            if nargin<2, B = 1; end
            T = A.cpop(B,'comp');
        end
        function B = sort(A)
           %SORT sort array of qtpl elements by frequecny
           w = [A.frequency];
           [~,I] = sort(w);
           B = A(I);
        end
        function [tpl,par] = get(obj,idx,w)
           %GET returns required tpl and par 
           
           if ~all(mod(idx,1)==0 & idx>0)
               error('second argument must be a vector of positive integers')
           end
           
           if nargin<3
               w=[]; 
           elseif ~(isscalar(w) & isnumeric(w))
               error('third input must be a numeric scalar')
           end
               
           if length(obj)>1 && isempty(w)
               error('for an array of template a third argument (frequency) is required');
           elseif length(obj)>1 
               Freq = [obj.frequency];
               T = obj(Freq==w);
               if isempty(T), error('requested frequency is not in template'); end
           else
               T = obj;
           end
           
           tpl = T.template(idx);
           par = T.parameters(:,idx);
           
        end
        function varargout = table(obj,idx,w)
            %DISP Displays the magnitude and paramerters in command window
            
            [tpl,par] = get(obj,idx,w);
            scases = sprintf('i_%i ',idx);
            ccases = split(scases(1:end-1)); % ignore 1 extra space at end
            
            if isempty(obj(1).parNames)
                spars = sprintf('par%i ',1:size(obj(1).parameters,1));
                cpars = split(spars(1:end-1)); % ignore 1 extra space at end
            else
                cpars = obj(1).parNames;
            end
            
            
            data = [real(tpl).' ; imag(tpl).' ; par];
            partab = array2table(data,'VariableNames',ccases,...
                'RowNames',{'deg', 'dB',cpars{:}});    
            if nargout==0
                disp(partab)
            elseif nargout==1
                varargout{1} = partab;
            else
                erorr('too many output argumetns')
            end
        end  
        function B = unwrap(A)
            %UNWRAP unwrap phase in qtpl
            %
            N = length(A);
            nom_qfr = nom(A);
            nomTpl = nom_qfr.response;
            unom = unwrap(real(nomTpl)*pi/180)*180/pi + 1i*imag(nomTpl);
            B  = A;
            for k=1:N
                nomPhase = real(unom(k));
                meanPhase = mean(real(A(k).template(2:end)));
                r = round((nomPhase - meanPhase)/360);
                if unom(k) == A(k).template(1)
                    B(k).template(2:end) = A(k).template(2:end)+r*360;
                else
                    B(k).template = A(k).template+r*360;
                end
            end
            
        end
        function G = nom(obj)
            %NOM returns a qfr object constracted of the nominal tpl points 
            N=length(obj);
            w = [obj.frequency].';
            g = zeros(N,1);
            for k=1:N
                g(k) = obj(k).template(1);
            end
            G = qfr(g,w);
        end
    end
    
    methods(Hidden = true)
        T = tplop( A,B,op );
    end
    
    methods(Static)
        %h = bodeplotter(tpl,w,opt,col)  % moved to utilities
        T = tplfile_import(filename,varargin);
    end
    
end

