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
            if nargin<4, options=[]; end
            
            switch method
                case 'grid', obj=obj.cgrid(w,0);
                case 'rndgrid', obj=obj.cgrid(w,1);
                case 'random', obj=obj.cgrid(w,2);
                case 'recgrid', obj=obj.adgrid(w,options);
                case 'aedgrid', obj=obj.adedge(w,options);
                otherwise, error('unrecognized method!')
            end
            
        end
        function obj = adgrid(obj,w,options)
            %ADGRID computes template by the recursive grid method
            
            % TODO: add  unstructured uncertainty 
            error('to do...')
            
        end
        function obj = adedge(obj,w,options)
            
            
            fprintf('Calculating templates by recurcive edge grid\n')
            
            Tacc=[5 2];
            
            npar=length(obj.pars);
            qdist = ([obj.pars.upper]' - [obj.pars.lower]')./double([obj.pars.cases]'); 
            indgrid = obj.idxgrid(npar);
            qg = grid(obj.pars,0,2);
            
            f = obj.qplant2func();

            %C = num2cell(qg,2);
            %C{end+1} = 0;
            c = qplant.pack(qg);
            
            Qpar = qg;
            tpl = qtpl(length(w)); % pre-allocating 
            col = distinguishable_colors(length(w));
            for kw = 1:length(w)
                fprintf('--> for w=%g [rad/s] \n',w(kw));
                s = 1j*w(kw);
                Tg = qplant.funcval(f,c,s);
                T = [];
                for k=1:npar
                    ind0=find(~indgrid(k,:));
                    ind1=find(indgrid(k,:));
                    for l=1:length(ind0)
                        Td=Tg(ind1(l))-Tg(ind0(l));
                        Td=rem(real(Td)+540,360)-180+1j*imag(Td);
                        qfix=max(qg(:,ind1(l))-qg(:,ind0(l)))/qdist(k);
                        if abs(real(Td)/Tacc(1)+1j*imag(Td)/Tacc(2))>=1 || qfix>1
                            [T1,Q1]=obj.recedge(f,s,qg(:,ind0(l)),qg(:,ind1(l)),Tg(ind0(l)),Tg(ind1(l)),Tacc,qfix);                           
                            T=[T T1];
                            Qpar=[Qpar Q1];        
                        end
                    end
                end
                scatter(real(T),imag(T),5,col(kw,:)); hold on
                tpl(kw) = qtpl(w(kw),T,Qpar);
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
            pck = obj.pack(pgrid);
            
            col = distinguishable_colors(length(w)); % to remove
            tpl = qtpl(length(w)); % pre-allocating 
            for k = 1:length(w)
                fprintf('--> for w=%g [rad/s] \n',w(k));
                %C{end}=1j*w(k);
                %nyq = f(C{:});               % in complex form
                T = obj.funcval(f,pck,1j*w(k));
                tpl(k)=qtpl(w(k),T,pgrid);  
                scatter(real(T),imag(T),10,col(k,:)); hold on
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
        function idx = idxgrid(N)
            %IDXGRID creates a grid of 0/1 
            %subroutine for ADEDGE
            j1=0:2^N-1;
            idx=zeros(N,2^N);
            for k=1:N
                idx(k,:)=rem(fix(j1/2^(k-1)),2);
            end
        end
        function c = pack(par)
            c = num2cell(par,2);
            c{end+1}=0;
        end
        function T = funcval(f,c,s)
            c{end} = s;
            nyq  = f(c{:});
            T=c2n(nyq,'unwrap');
        end
        function [Tnew,Qpar] = recedge(trf,s,qmin,qmax,Tmin,Tmax,Tacc,qdist)
            %RECEDGE    subroutine used by ADEDGE
            %           [Tnew,Qpar]=recedge(trf,s,qmin,qmax,Tmin,Tmax,Tacc,qdist);
            %
            %   Adopted from Qsyn: Original Author: M Nordin
            
            %global prune_on
            qbreak=abs(qmin-qmax);
            qbreak=min(qbreak((qbreak~=0)));
            if qbreak<1e-8
                prune_on=0; % return value???
                return
            end
            qnew=0.5*(qmin+qmax);
            qdist=0.5*qdist;
            
            %avoid NaN from ~plant.m = xxxplant.m  peo 970318
            %value = feval(trf,qnew,s);
            c=qplant.pack(qnew);
            Tcurr = qplant.funcval(trf,c,s);
            if isnan(Tcurr) 
                %value = feval(trf,qnew+eps*ones(size(qnew)),s); 
                c=obj.pack(qnew+eps*ones(size(qnew)));
                Tcurr = obj.funcval(trf,c,s);
            end
            
            % Tcurr=c2n(feval(trf,qnew,s));
            %Tcurr=c2n(value,-180);  %peo960705
            Qpar=qnew;
            Tnew=Tcurr;
            Td=Tcurr-Tmin;
            Td=rem(real(Td)+540,360)-180+1j*imag(Td);
            %The wrapped distance
            if abs(real(Td)/Tacc(1)+1j*imag(Td)/Tacc(2))>1 || qdist>=1
                [T1,Q1]=qplant.recedge(trf,s,qmin,qnew,Tmin,Tcurr,Tacc,qdist);
                Tnew=[Tnew T1];
                Qpar=[Qpar Q1];
            end
            Td=Tmax-Tcurr;
            Td=rem(real(Td)+540,360)-180+1j*imag(Td);
            if abs(real(Td)/Tacc(1)+1j*imag(Td)/Tacc(2))>1 || qdist>=1
                [T1,Q1]=qplant.recedge(trf,s,qnew,qmax,Tcurr,Tmax,Tacc,qdist);
                Tnew=[T1 Tnew];
                Qpar=[Q1 Qpar];
            end
        end
    end
        
end

