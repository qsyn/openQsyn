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
    end
    
end

