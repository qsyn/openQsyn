classdef qmeas < qplant
    %QMEAS measurement based qplant 
    %   
    % Serves to perform bound calculations when only frequency response
    % data is aviable
    
    properties
        
    end
    
    methods
        function obj = qmeas(tpl)
            %QMEAS construct an instance of this class
            %
            %   P = QMEAS(T)   creates a qmeas object P from the qtpl object T
            %
          
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

