classdef qplant < handle
    %QPLANT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        num
        den
        pars
        templates
        nom
        info
    end
    
    methods
        function obj = qplant(num,den)
            %QPLANT Construct an instance of this class
            %   Detailed explanation goes here
            obj.num = num;
            obj.den = den;
            npar=[];
            dpar=[];
            if isa(num,'qpoly'), npar = num.pars; end
            if isa(den,'qpoly'), dpar = den.pars; end
            obj.pars = unique(vertcat(npar,dpar));
            
            obj.info=sprintf('generated from [num,den] data on: %s',datetime);
        end
        function obj = cnom(obj,w)
            %CPNOM computes the nominal transfer function
            %   Detailed explanation goes here
           if nargin<2, w = logspace(-1,2,50); end 
           f = qplant2func(obj);
           pnom = [obj.pars.nominal];
           C = num2cell(pnom);
           C{end+1}=1j*w;
           nyq = f(C{:});               % in complex form
           obj.nom=qfr(c2n(nyq),w) ;
        end  
        function obj = ctpl(obj,method,w,options)
            %CTPL computes the templates...
            
            
        end
        function obj = cgrid(obj,w)
            %CGRID computes tpl by the grid method
            
            % TODO: add  unstructured uncertainty 
            
            disp('Calculating templates using the grid method ');
            f = qplant2func(obj);
            pgrid = grid(obj.pars);
            C = num2cell(pgrid,2);
            C{end+1}=0;
            col = distinguishable_colors(length(w));
            for k = 1:length(w)
                fprintf('--> for w=%g [rad/s] \n',w(k));
                C{end}=1j*w(k);
                nyq = f(C{:});               % in complex form
                tpl(k)=qtpl(w,c2n(nyq),pgrid);  
                scatter(real(c2n(nyq)),imag(c2n(nyq)),10,col(k,:)); hold on
            end 
            obj.templates = tpl;
        end
        function h = qplant2func(obj)
            snum = obj.poly2str(obj.num);
            sden = obj.poly2str(obj.den);
            
            args = sprintf('%s, ',obj.pars.name);
            argF = sprintf('@(%s, s) ',args(1:end-2));
            s = sprintf('%s (%s)./(%s)',argF,snum,sden);
         
            h = str2func(s);
        end
        function P = tf(obj)
            %TF convert Qplant to its NOMINAL transfer function.
            
            if isnumeric(obj.num)
                NUM = obj.num;
            else
                NUM = obj.num.nom;
            end
            if isnumeric(obj.den)
                DEN = obj.den;
            else
                DEN = obj.den.nom;
            end
            P = tf(NUM,DEN);
        end
    end
    
    methods(Static)
        function s = poly2str(p)
            if isa(p,'qpoly')
                fnum = qpoly2func(p);
                snum = func2str(fnum);
                inum = strfind(snum,')');
                s = snum(inum+1:end);
            else % is polynom
                s ='';
                for k=length(p):-1:1
                    s = sprintf('%s %g*s.^%i + ',s,p(k),k-1);
                end
                s = s(1:end-2);
            end
            
        end
        


                
    end
        
end

