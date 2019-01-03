classdef qdesign  < handle
    %QLD is a class for QFT Loop Design -- A class that host all methods
    %required for SISO QFT loop shaping design
    
    properties
        tpl
        nom
        spc
        bnd
    end
    
    methods
        function obj = qdesign(plant,specs)
            %QLD Construct an instance of the QLD class
            %   Detailed explanation goes here
            obj.tpl = plant.templates;
            obj.nom = plant.nominal;
            obj.spc = specs;
            disp(['You now have a QFT loop desgin object. ',... 
                 'Compute bounds using CBND'])
        end
    end
    
    methods(Static)
        bound = makebnd(tpl,specfunc,spec,gphase,gmag)
    end
end

