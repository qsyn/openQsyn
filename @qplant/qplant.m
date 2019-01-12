classdef qplant < handle
    %QPLANT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        num
        den
        pars
        templates
        nominal
        info
    end
    
    methods
        function obj = qplant( num,den ) 
            %QPLANT Construct an instance of the QPLANT class
            %   Detailed explanation goes here
            %
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
           nyq = f(C{:});                   
           obj.nominal=qfr(c2n(nyq,'unwarp'),w) ;
        end  
        function obj = ctpl(obj,method,w,options)
            %CTPL computes the templates for given qplant object
            %
            %   obj = ctpl(obj,method,w,options)
            %
            %   Inputs
            %
            %   obj     qplant object
            %
            %   method  template computation method. one of the following:
            %               'grid'      uniform grid over parameter
            %               'rngrid'    random grid over the parameters
            %               'random'    random sample points 
            %               'recgrid'   recurcive grid (TO DO)
            %               'aedgrid'   recurcive edge grid
            % 
            %   w       frequency vector
            %
            %   options CURRENTLY NO USAGE
            
            % TODO: add additional options...
            if nargin<4, options=[]; end
                        
            switch method
                case 'grid', tpl=obj.cgrid(w,0);
                case 'rndgrid', tpl=obj.cgrid(w,1);
                case 'random', tpl=obj.cgrid(w,2);
                case 'recgrid', tpl=obj.adgrid(w,options);
                case 'aedgrid', tpl=obj.adedge(w,options);
                otherwise, error('unrecognized method!')
            end

            % add nominal point at beginning of each template
            pnom = [obj.pars.nominal]';
            tnom=cases(obj,pnom,w);
            for k=1:length(w)
                tpl(k)=add2tpl(tpl(k),tnom(k),pnom(:,1),'x');
            end
            
            if isempty(obj.templates)
                obj.templates = tpl;
            else 
                % replacing existing freqeucies and inserting new ones in
                % the right places (sorted by frequency)
                w0 = [obj.templates.frequency];
                w1 = [tpl.frequency];
                w = unique([w0 w1]);
                inew = ismember(w,w1);  
                i0 = ismember(w(~inew),w0);
                TPL = qtpl(length(w));
                TPL(inew) = tpl;
                TPL(~inew) = obj.templates(i0);
                obj.templates = TPL;
            end
          
            
        end
        function tpl = adgrid(obj,w,options)
            %ADGRID computes template by the recursive grid method
            
            % TODO: add  unstructured uncertainty 
            tpl = [];
            error('to do...')
            
            
        end
        function tpl = adedge(obj,w,options)
            
            
            fprintf('Calculating templates by recurcive edge grid\n')
            
            plot_on = 0; % to be given as option in future
            Tacc=[5 2];  % to be given as option in future
            
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
              
                idx=boundary(real(T)',imag(T)',0.3); % replaces PRUNE (introduced in R2014b)
                if plot_on
                    %scatter(real(T),imag(T),5,col(kw,:)); hold on
                    %scatter(real(T(idx)),imag(T(idx)),10,col(kw,:),'marker','o');
                end
                tpl(kw) = qtpl(w(kw),[T(idx).'],Qpar(:,idx));
            end
            
        end
        function tpl = cgrid(obj,w,rnd)
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
                tpl(k)=qtpl(w(k),T.',pgrid);  
                scatter(real(T),imag(T),10,col(k,:)); hold on
            end 
        end
        function h = qplant2func(obj)
            %QPLANT2FUNC return an handle to a function object f@(p1,p2,...pn,s) 
            %with p1,...,pn corresponding to uncertain parameters.
            %
            %Used by: FUNCVAL, CASES, CTPL, CNOM and their subroutines
            
            snum = obj.poly2str(obj.num);
            sden = obj.poly2str(obj.den);
            
            args = sprintf('%s, ',obj.pars.name);
            argF = sprintf('@(%s, s) ',args(1:end-2));
            s = sprintf('%s (%s)./(%s)',argF,snum,sden);
         
            h = str2func(s);
        end
        function P = tf(obj,par)
            %TF converts Qplant to a transfer function
            %
            %   P = TF(QPLANT)   returns the transfer function for the 
            %   nominal case
            % 
            %   P = TF(QPLANT,PAR)   returns the transfer function for the
            %   case given by vector of parameters PAR
            if nargin<2, par = []; end
            if ~isempty(par) && length(par)~=length(obj.pars)
                error('Number of input parameters must be same as in qplant object')
            end
            
            if isnumeric(obj.num)
                NUM = obj.num;
            elseif isempty(par)
                NUM = obj.num.nom;
            else
                ip = ismember(obj.pars,obj.num.pars);
                npar = par(ip);
                NUM = obj.num.polycase(npar);
            end
            if isnumeric(obj.den)
                DEN = obj.den;
            elseif isempty(par)
                DEN = obj.den.nom;
            else
                ip = ismember(obj.pars,obj.den.pars);
                dpar = par(ip);
                DEN = obj.den.polycase(dpar);
            end
            P = tf(NUM,DEN);
        end
        function varargout = showtpl(obj,w,varargin)
            %SHOWTPL plots the templates at given frequencies 
            %This is basically a wrapper for qtpl.show to allow more
            %conviniante access and allow the original Qsyn capabilities
            %
            %  	showtpl(QPLANT)     displays template QTPL 
            %   
            %
            %   showtpl(QPLANT,W)   display only tempaltes at freqeuncies W       
            %   
            %   showtpl(TPLF,W,FHAND)    draws in figure with handle FHAND   
            %
            %   showtpl(TPLF,W,PARAMETER,VALUE)   use parameter/value pairs to
            %   specify additional properties:
            %       'mode'    string specifing the mode
            %       'color'   color array in RGB format
            %       'marker'  string for marker points
            %      	'fill' 	  boolian scalar 1 | 0 (def)
            %       'case'    vector of indices specifying case(s) to show 
            %
            %   QPLANT.showtpl(...)     alternative usage
            %
            %   Avilable modes:
            %       'nom'(def)	The nominal plant is displayed and the
            %                   templates are drawn correctly relative 
            %               	their nominal points 
            %       'point' 	The user clicks with his mouse on the Nichols
            %                  	chart for the location of the nominal point 
            %                   of the next template --> NOT IMPLEMENTED!!!
            %       'nonom'     plot templates in their position w/o nominal
            
            %%% Input handling
            if nargin<2, w = []; end           
            
            wtpl = [obj.templates.frequency];
            if isempty(w), w = wtpl; end
                       
            ishow = ismember(wtpl,w);
            if all(~ishow)
                error('w must be a subset of the avialble frequencies'); 
            end              
            opt = 'nom';
            if nargin>2
                k=0;
                if ishandle(varargin{1}), k=1; end
                imode = strcmp({varargin{(1+k):2:end}},'mode');
                if any(imode)
                     opt = varargin{2*find(imode)+k};
                     idx = true(1,length(varargin));
                     idx(k+[2*find(imode)-1:2*find(imode)]) = false(1,2);
                     ARG = {varargin{idx}};
                else
                    ARG = varargin;
                end
            end
            
            %%% Main 
            h = obj.templates(ishow).show(ARG{:});
            
            % plot nominal
            if strcmp(opt,'nom')
                if isempty(obj.nominal)
                    disp('no nominal exists.');
                else
                    %plot(real(obj.nominal),imag(obj.nominal));
                    obj.nominal.show();
                end
            elseif ~strcmp(opt,'nonom')
                error('mode options: ''nom'' | ''nonom''')
            end
           
            if nargout==1
                varargout{1}=h;
            end
        end
        function bodcases(obj,par,w,opt)
            %BODCASES plant frequency domain simulation for user selected
            %cases on Bode plont
            %
            %   bodcases(QPLANT)   plots bode for all cases given by the 
            %   plant parameters
            %
            %   bodcases(QPLANT,W)   specify the frequecnies
            %
            %   bodcases(QPLANT,W,OPT)   specify wihch part to plot: 
            %                            'mag' | 'phase' | 'magphase' (def)
            %
            %   
            %   BODCASES and NICCASES replace CASES in Qsyn
            
            %   TODO: something with OPT
            if nargin<4, opt=[]; end                    
            if nargin<3, w=[]; end
            if nargin<2, par=[]; end
            
            if isempty(opt), opt = 'magphase'; end    
            
            [res,w] = cases(obj,par,w);
            
            col = distinguishable_colors(size(res,2));
            qtpl.bodeplotter(res.',w,opt,col);           
                
        end
        function varargout=cases(obj,par,w)
            %CASES returns the template points for given parametric cases
            %It does not plot anything!
            
            if nargin<3, w=[]; end
            if nargin<2, par=[]; end
            
            if isempty(w), w = obj.nominal.frequency; end
            w = reshape(w,[],1); % make sure w is a column vector.
            
            f = obj.qplant2func();
            if isempty(par)
                pgrid = grid(obj.pars,0);
            else
                pgrid = par;
            end
            p = qplant.pack(pgrid);
            
            varargout{1} = obj.funcval(f,p,w*1j); % response in Nichols format
            if nargout==2
                varargout{2} = w;
            elseif nargout>2
                error('to many outputs!')
            end
        end
        function step(obj,pars,tfinal)
            %STEP computes step response for plant cases
            if nargin<3, tfinal=[]; end
            if nargin<2, pars=[]; end
            
            if isempty(pars)
                pgrid = grid(obj.pars,0);
            else
                pgrid = pars;
            end

            for k=1:size(pgrid,2)
                Pk = tf(obj,pgrid(:,k));
                [y,t]=step(Pk,tfinal);
                plot(t,y);
                hold on
            end
            xlabel('Time [s]');
            ylabel('Amplitude');
            axis tight
            hold off           
        end
        function [tpl,par] = gettpl(obj,idx,w)
            %GETTPL get tpl point(s) from a plant 
            Freq = [obj.templates.frequency];
            T = obj.templates(Freq==w);
            if isempty(T), error('requested frequency is not in template'); end
            
            [tpl,par] = get(T,idx);
            
        end
    end
    
    methods(Static)
        function s = poly2str(p)
            %POLY2STR converts a polynom to a string
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
            %PACK returns a cell array containning all parameters in PAR 
            %and leaves one additional spot empty
            c = num2cell(par,2);
            c{end+1}=0;
        end
        function T = funcval(f,c,s)
            %FUNCVAL return the value of the plant, represented by a 
            %function handle F created by QPLANT2FUNC, for given parameter 
            %case C and for S. The output is in unwrapped Nihols format.
            c{end} = s;
            %c = num2cell(par,2);
            %c{end+1} = s;
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

