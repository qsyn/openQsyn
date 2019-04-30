classdef qmeas < qplant
    %QMEAS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = qmeas(tpl)
            %QBLACKBOX Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@qplant(tpl);
        end
        
        function cnom(varargin)
            %CNOM computes nominal case 
            disp('cannot perform on measured plants. use INTERPNOM instead.')
        end
        function ctpl(varargin)
            disp('meathod not avilable for measured plants')
        end
    end
    
end

