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
            
            obj.frequency = double(frequency);
            obj.template = double(template);
            obj.parameters = double(parameters); 
        end  
    end
    
end

