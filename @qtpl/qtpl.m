classdef qtpl
    %QFTTPL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frequency
        template
        parameters
    end
    
    methods
        function obj = qtpl(frequency,template,parameters)
            %QTPL construct a qtpl object
            if nargin==0, return; end
            if nargin==1 % build an empty array of length given by input 1
                obj(frequency,1) = obj;
                return
            end
            obj.frequency = double(frequency);
            obj.template = double(template);
            obj.parameters = double(parameters); 
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
           par = T.parameters(idx);
           
        end
    end
    
    methods(Static)
        h = bodeplotter(tpl,w,opt,col) 
        T = tplfile_import(filename)
    end
    
end

