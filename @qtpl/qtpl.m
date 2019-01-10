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
    end
    
    methods(Static)
        h = bodeplotter(tpl,w,opt,col) 
    end
    
end

