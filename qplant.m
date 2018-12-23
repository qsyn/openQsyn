classdef qplant
    %QPLANT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Num
        Den
        Parameters
        Templates
        Nominal
        Info
    end
    
    methods
        function obj = qplant(num,den)
            %QPLANT Construct an instance of this class
            %   Detailed explanation goes here
            obj.Num = num;
            obj.Den = den;
        end
        
        function obj = gpars(obj)
            % get parameters from num and den
            
            pars
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

