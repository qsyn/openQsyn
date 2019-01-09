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
        function h = loopnic(obj,C)
            %LOOPNIC plots the open-loop on a Nichols chart
            
            % defaults
            linecolor = 'k'; 
            t_color = distinguishable_colors(length(obj.tpl)); 
            hold on
            
            Lnom = series(obj.nom,C);
            h = Lnom.show('color',linecolor);
            
            for k=1:length(obj.tpl)
                tk = qfr(obj.tpl(k).template(1),obj.tpl(k).frequency);
                Ltpl = series(tk,C);
                h(end+1) = show(Ltpl,'marker','square',...
                     'markeredgecolor','k','markerfacecolor',t_color(k,:));
            end
            
        end
    end
    
    methods(Static)
        bound = makebnd(tpl,specfunc,spec,gphase,gmag)
    end
end

