classdef qdesign  < handle
    %QDESIGN is a class for QFT Loop Design -- A class that hosts all 
    %methods required for SISO QFT loop shaping design
    
    properties
        tpl         % templates (qtpl array)
        nom         % nominal case (qfr)
        spc         % specifications (qspc array)
        bnd         % bounds (struct)
        col         % color for each frequency (rgb array)
    end
    
    methods
        function obj = qdesign(plant,specs)
            %QLD Construct an instance of the QDESIGN class
            %   
            %Usage: 
            %
            %    obj = QDESIGN(plant,specs)    constract a QDESIGN object
            %    from given qplant and qspc
            %
            %Part of the Open Qsyn toolbox. 
            
            obj.tpl = plant.templates;
            obj.nom = plant.nominal;
            obj.spc = specs;
            obj.col = lines(length(obj.tpl)); % preserves color for each freq
            disp(['You now have a QFT loop desgin object. ',... 
                 'Compute bounds using CBND'])
        end
        function [] = clmag(obj,C,F)
            %CLMAG computes closed loop magnitude response
            %
            %   obj     qdesign object 
            %
            %   C       feedback compensator in LTI foramt or an abolute
            %           gain; Default C=1.
            %
            %   F       prefilter in in LTI foramt or an abolute
            %           gain. Default F=1.
            %
            %   Replaces Qsyn command FDESIGN
            
            if nargin<3, F=[]; end
            if isempty(F), F=1; end
            
            if isnumeric(C), C=tf(C); end
            if isnumeric(F), F=tf(F); end
            
            % compute nominal
            wnom = obj.nom.frequency;
            Cwnom = squeeze(freqresp(C,wnom)); % in complex plain
            Fwnom = squeeze(freqresp(F,wnom)); % in complex plain
            Pnom = n2c(obj.nom.response);
            Lnom = Cwnom.*Pnom;
            Tnom = Fwnom.*Lnom./(1+Lnom);
            
            qfr_nom = qfr(c2n(Tnom),wnom);
            bodeplot(qfr_nom,'PhaseVisible','off')
            hold on
            
            % compute templates
            Tu = cpop( obj.tpl,C,'comp');
            Tr = cpop( Tu,F,'*');
            
            [mag,~,w] = Tr.bode([]);
            mag_min = min(mag');
            mag_max = max(mag');
            scatter(w,mag_min,'bo')
            scatter(w,mag_max,'bo')
        end
        function varargin = loopnic(obj,C)
            %LOOPNIC plots the open-loop on a Nichols chart
            
            % defaults
            linecolor = 'k'; 
            %t_color = distinguishable_colors(length(obj.tpl)); 
            t_color = obj.col;
            hold on
            
            Lnom = series(obj.nom,C);
            h = Lnom.show('color',linecolor);
            
            for k=1:length(obj.tpl)
                tk = qfr(obj.tpl(k).template(1),obj.tpl(k).frequency);
                Ltpl = series(tk,C);
                h(end+1) = show(Ltpl,'marker','square',...
                     'markeredgecolor','k','markerfacecolor',t_color(k,:));
                text(real(Ltpl.response),imag(Ltpl.response),...
                    sprintf(' %g',obj.tpl(k).frequency),'clipping','on') % single space added
            end
            
            if nargout==1, varargin{1}=h; end
            
        end
    end
    
    methods(Static)
        bound = makebnd(tpl,specfunc,spec,gphase,gmag)
    end
end

