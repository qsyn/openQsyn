classdef qtpl
    %QFTTPL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frequency
        template
        parameters
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
            
        end  
        function [ T ] = plus( A,B )
            %PLUS adds two qtpl arrays
            %
            %   [ T ] = plus( A,B )   performs an addition between qtpl objects A and
            %   B, returns output as a qtpl object T.
            %
            %   Note that the plus opperation is performed in Nichols form (deg+i*db),
            %   i.e., for siso transfer functions A,B: plus(A,B) = A(s)*B(s).
            T = tplop(A,B,'+');
        end
        function [ T ] = minus( A,B )
            %PLUS substruct two qtpl arrays
            %
            %   [ T ] = minus( A,B )   performs a substraction between qtpl
            %   objects A and B, returns output as a qtpl object T.
            %
            %   Note that the plus opperation is performed in Nichols form (deg+i*db),
            %   i.e., for siso transfer functions A,B: minu(A,B) = A(s)/B(s).
            T = tplop(A,B,'-');
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
        function table(obj,idx,w)
            %DISP Displays the magnitude and paramerters in command window
            
            [tpl,par] = get(obj,idx,w);
            scases = sprintf('i_%i ',idx);
            ccases = split(scases(1:end-1)); % ignore 1 extra space at end
            
            spars = sprintf('par%i ',1:size(obj(1).parameters,1));
            cpars = split(spars(1:end-1)); % ignore 1 extra space at end
            
            data = [real(tpl).' ; imag(tpl).' ; par];
            partab = array2table(data,'VariableNames',ccases,...
                'RowNames',{'deg', 'dB',cpars{:}});    

            disp(partab)
        end
        function obj = unwrap(obj)
            %UNWRAP unwrap phase in qtpl
            %
            N = length(obj);
            nom = get(obj,1);
            unom = unwrap(nom*pi/180)*180/pi + 1i*imag(nom);
            for k=1:N
                nomPhase = real(unom(k));
                meanPhase = mean(real(obj(k).template(2:end)));
                r = round((nomPhase - meanPhase)/360);
                obj(k).template(2:end) = obj(k).template(2:end)+r*360;
            end
            
        end
    end
    
    methods(Static)
        h = bodeplotter(tpl,w,opt,col) 
        T = tplfile_import(filename,varargin)
    end
    
end

