classdef qfr
    %QFR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frequency
        nic
    end
    
    methods
        function obj = qfr(nresponse,frequecny)
            %QFR Construct an instance of this class
            %   Detailed explanation goes here
            obj.frequency = frequecny;
            obj.nic = nresponse;
        end
        function [] = show(obj)
            %NIC plots the on Nichols chart
            plot(real(obj.nic),imag(obj.nic),'-k','linewidth',1.5)
        end
        function [] = nichols(obj)
            %NICHOLS plots a nichols chart, compatible to Control System
            %Toolbox
            G = obj.frd();
            nichols(G,obj.frequency);
        end
        function [] = bode(obj)
            %BODE plots a bode plot
            G = obj.frd();
            bode(G,obj.frequency);
        end
        function G = frd(obj)
            %FRD converts qfr to frd object
            
            response = n2c(obj.nic);
            G = frd(response,obj.frequency);    
        end 
    end
end

