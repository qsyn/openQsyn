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
            
            % TODO: add additional options...
            switch method
                case 'grid', obj=obj.cgrid(w,0);
                case 'rndgrid', obj=obj.cgrid(w,1);
                case 'random', obj=obj.cgrid(w,2);
                otherwise, error('unrecognized method!')
            end
            
            
        end
        function obj = cgrid(obj,w,rnd)
            %CGRID computes tpl by the grid method
            %   facilitates grid, random grid, and random sampling 
            %   
            %   for random and random grid the parameter set is random, yet
            %   identical in every frequency.
            
            % TODO: add  unstructured uncertainty 
            %       add  options to change cases 
            if nargin<3, rnd=0; end
            switch rnd
                case 0
                    method='grid'; 
                    pgrid = grid(obj.pars,rnd);
                case 1
                    method='random grid'; 
                    pgrid = grid(obj.pars,rnd);
                case 2 
                    method='random sampling';
                    pgrid = sample(obj.pars,100); % correct usage: options.cases(=100)
                otherwise, error('unavilable rnd option')
            end
                
            fprintf('Calculating templates using the %s method \n',method)
            f = qplant2func(obj);
            
            C = num2cell(pgrid,2);
            C{end+1}=0;
            col = distinguishable_colors(length(w)); % to remove
            tpl = qtpl(length(w)); % pre-allocating 
            for k = 1:length(w)
                fprintf('--> for w=%g [rad/s] \n',w(k));
                C{end}=1j*w(k);
                nyq = f(C{:});               % in complex form
                tpl(k)=qtpl(w(k),c2n(nyq),pgrid);  
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

