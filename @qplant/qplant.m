classdef qplant < handle
    %QPLANT class for plant description and related operations
    %
    %   This class contains the plant description and all methods related to
    %   the plant:
    %   - freqeuncy response for uncertain cases
    %   - time domain simulations for uncertain cases
    %   - extraction of parameters and properties
    %   - template computation
    %
    %Construction:
    %   P = QPLANT(num,den)     constructs a qplant object P from numerator and
    %   denominator data - both are qpoly objects.
    %
    %Adding a delay:
    %   P.adelay(h)     adds a delay h to the plant. the delay may be a constant
    %   scalar or uncertain element (qpar or qexpression)
    %
    %Adding unstructured uncertainty:
    %   P.aunstruc(w,m) adds an unstructured uncertainty given by a vector of
    %   absolute magnitudes corrsponding to given frequency vectors
    %
    %Plotting a Bode plot:
    %   P.bodecases(w,par)   plots a Bode plot for frequencies in w and
    %   parameter cases given in the matrix par.
    %
    %Computing templates
    %   P.ctpl(method,w)    computes tmeplate according to specified method at
    %   given frequencies w
    %
    %For a list QPLANT methods type: methods qplant
    %For help on a specific methods type: help qplant/<method>
    
    properties (SetAccess = 'protected')
        num             % num  (qpoly array)
        den             % den (qpoly array)
        pars            % uncertain parameters (qpar array)
        delay           % delay (string)
        unstruct        % unstructed uncertainty
        uncint          % uuncertain integrator/diffrentiator
    end
    
    properties
        info        char   % string
        templates   qtpl   % template (qtpl array)
        nominal     qfr    % nominal case (qfr)
    end
    
    
    methods
        function obj = qplant( varargin )
            %QPLANT Construct an instance of the QPLANT class
            %   P = QPLANT(num,den)     constructs a qplant object P from 
            %   numerator and denominator data
            %
            %   P = QPLANT(func_handle,pars)   constructs a 'black box' qplant
            %   object P from a function handle and a qpar array. 
            %
            %   P = QPLANT(tpl)     constructs a measurment based qplant
            %   from qtpl array

            % measured plant:
            if nargin==1
                if ~isa(varargin{1},'qtpl')
                    error('When called with a single input it has to be of qtpl class');
                end
                obj.templates = varargin{1};
                obj.info = 'Plat constracted from qtpl data';
                return
            end
            
            % num/den or blackbox plant:
            isBlackBox = 0; % default type
            % check first input
            switch class(varargin{1})
                case {'qpoly', 'qexpression'}
                    obj.num = varargin{1};
                    npar = varargin{1}.pars;
                case 'qpar'
                    obj.num = varargin{1};
                    npar = varargin{1};
                case {'double','single','int32','int64'}
                    validateattributes(varargin{1},{'numeric'},...
                        {'row','nonempty','real'},'qplant','num')
                    obj.num = varargin{1};
                    npar = [];
                case 'function_handle'
                    isBlackBox = 1;
                    obj.blackBox = varargin{1};
                    npar =[];
                otherwise
                    error('input 1 is not a supported class.')
            end
            % check second input
            switch class(varargin{2})
                case {'qpoly', 'qexpression'}
                    if isBlackBox
                        error('for a black box model the 2nd argument must be a qpar array');
                    end
                    obj.den = varargin{2};
                    dpar = varargin{2}.pars;
                case 'qpar'
                    if ~isBlackBox
                        obj.den = varargin{2};
                    end
                    dpar = varargin{2};
                case {'double','single','int8','int16','int32','int64'}
                    if isBlackBox
                        error('for a black box model the 2nd argument must be a qpar array');
                    end
                    validateattributes(varargin{2},{'numeric'},...
                        {'row','nonempty','real'},'qplant','num')
                    obj.den = varargin{2};
                    dpar = [];
                otherwise
                    error('input 2 is not a supported class.')
            end
            
            % build list of qpars
            obj.pars = unique(vertcat(npar,dpar));
            
            %obj.info=sprintf('generated from [num,den] data on: %s',datetime);%cannot work on older version
            if isBlackBox
                obj.info=sprintf('generated from BalckBox func handle on: %s',datestr(datetime));
            else
                obj.info=sprintf('generated from [num,den] data on: %s',datestr(datetime));
            end
            
        end
        function obj = adelay(obj,del)
            %ADELAY adds a delay to an existing plant without delay
            if ~isempty(obj.delay), error('plant already have a delay'); end
            if isa(del,'qexpression')
                obj.delay = del;
                obj.pars = unique(vertcat(obj.pars ,del.pars));
            elseif isa(del,'qpar')
                obj.delay = del;
                obj.pars = unique(vertcat(obj.pars ,del));
            elseif isnumeric(num) && isscalar(num)
                obj.delay = del;
            else
                error('unsupported object for delay')
            end
            
            
        end
        function obj = auncint(obj,pow)
            %Auncint adds uncertain integrators/diffrentiators
            
            if ~(isnumeric(pow) && isvector(pow) && all(mod(pow,1)==0))
                error('second input argument must be a vector of integers')
            end
            pow0=unique([0 pow]);
            uncint_par = qpar('uncint_par',0,pow0,'disc');
            
            if isempty(obj.uncint)
                obj.pars = vertcat(obj.pars, uncint_par);
            else
                idx = strcmp({obj.pars.name},'uncint_par');
                obj.pars(idx) = uncint_par;
            end
            obj.uncint = pow;
            
        end
        function obj = aunstruc(obj,w,absval)
            %AUNSTRUCT addds plant unsturcured uncertainty
            p= inputParser;
            addRequired(p,'w',@(x)validateattributes(x,{'numeric'},{'nonempty','vector','positive','real'}));
            addRequired(p,'absval',@(x)validateattributes(x,{'numeric'},{'nonempty','vector','nonnegative'}));
            parse(p,w,absval);
            w = reshape(p.Results.w,[],1);
            absval = reshape(p.Results.absval,[],1);
            obj.unstruct = [w absval];
        end
        function obj = cnom(obj,w)
            %CPNOM computes the nominal transfer function
            %   
            %   CNOM(P,W)   computes the nominal for qplant object P for the frequency 
            %   vector W; results is stored under the property 'nominal'.
            %
            %   See also: qplant/ctpl
            if nargin<2, w = logspace(-1,2,50); end
            f = qplant2func(obj);
            pnom = [obj.pars.nominal];
            C = num2cell(pnom);
            C{end+1}=1j*w;
            nyq = f(C{:});
            obj.nominal=qfr(c2n(nyq,'unwarp'),w) ;
        end
        function tpl = cgrid(obj,w,rnd,n)
            %CGRID computes tpl by the grid method
            %   facilitates grid, random grid, and random sampling
            %
            %   for random and random grid the parameter set is random, yet
            %   identical in every frequency.
            %
            if nargin<3, rnd=0; end
            idx = ~strcmp({obj.pars.name},'uncint_par');
            switch rnd
                case 0
                    method='grid';
                    pgrid = grid(obj.pars(idx),n,rnd);
                case 1
                    method='random grid';
                    pgrid = grid(obj.pars(idx),n,rnd);
                case 2
                    method='random sampling';
                    pgrid = sample(obj.pars(idx),n); % correct usage: options.cases(=100)
                otherwise, error('unavilable rnd option')
            end
            
            fprintf('Calculating templates using the %s method \n',method)
            f = qplant2func(obj);
            pck = obj.pack(pgrid);
            
            %col = distinguishable_colors(length(w)); % to remove
            %col = lines(length(w));
            tpl = qtpl(length(w)); % pre-allocating
            for k = 1:length(w)
                fprintf('--> for w=%g [rad/s] \n',w(k));
                %C{end}=1j*w(k);
                %nyq = f(C{:});               % in complex form
                T = obj.funcval(f,pck,1j*w(k));
                tpl(k)=qtpl(w(k),T.',pgrid);
                %scatter(real(T),imag(T),10,col(k,:)); hold on
            end
        end
        function tpl = addunstruct(obj,tpl_in)
            %CUNSTRUCT adds unstructured uncertainty into template
            w = obj.unstruct(:,1);
            absmag = obj.unstruct(:,2); % magintude (abs)
            tpl = qtpl(length(w));
            a = linspace(0,2*pi,7).';
            a = a(1:end-1);
            for k=1:length(w)
                if absmag(k) > 0
                    circ = absmag(k)*(cos(a) + 1j*sin(sin(a))); % cirlce in complex plain
                    tpl(k) = cpop(tpl_in(k),circ,'+');
                else
                    tpl(k) = tpl_in(k);
                end
            end
            
        end
        function tpl = cases2tpl(obj,par,w)
            %CASES2TPL creates templates based on given by parameter cases
            %
            if size(par,1) ~= length(obj.pars)
                error('par must have %i raws', length(obj.pars));
            end
            p = obj.pack(par);
            f = obj.qplant2func();
            tpl = qtpl(length(w)); % pre-allocating
            for k=1:length(w)
                T = obj.funcval(f,p,1j*w(k));
                tpl(k)=qtpl(w(k),T.',par);
            end
        end
        function h = qplant2func(obj)
            %QPLANT2FUNC return an handle to a function object f@(p1,p2,...pn,s)
            %with p1,...,pn corresponding to uncertain parameters.
            %
            %Used by: FUNCVAL, CASES, CTPL, CNOM and their subroutines
            
            %if ~isempty(obj.blackBox)
            %    h = obj.blackBox;
            %    return
            %end
            % if not a black box continue as follows:
            
            p = obj.num;
            if isa(p,'qpoly') || ( isnumeric(p) && isrow(p) )
                snum = obj.poly2str(p);
            elseif isa(p,'qexpression')
                snum = strrep(strrep(strrep(p.expression, '*', '.*'),'/','./'),'^','.^');
            elseif isa(p,'qpar')
                snum = p.name;
            end
            
            p = obj.den;
            if isa(p,'qpoly') || ( isnumeric(p) && isrow(p) )
                sden = obj.poly2str(p);
            elseif isa(p,'qexpression')
                sden = strrep(strrep(strrep(p.expression, '*', '.*'),'/','./'),'^','.^');
            elseif isa(p,'qpar')
                sden = p.name;
            end
            
            if isempty(obj.delay)
                sdel = '';
            elseif isa(obj.delay,'qexpression')
                sdel = sprintf('.*exp(-s.*%s)',obj.delay.expression);
            elseif isa(obj.delay,'qpar')
                sdel = sprintf('.*exp(-s.*%s)',obj.delay.name);
            elseif isnumeric(obj.delay) && isscalar(obj.delay)
                sdel = sprintf('.*exp(-s.*%g)',obj.delay);
            else
                error('what kind of a delay is that?!')
            end
            
            if isempty(obj.uncint)
                sint = '';
            else
                sint = sprintf('(s.^uncint_par).*');
            end
            
            args = sprintf('%s, ',obj.pars.name);
            argF = sprintf('@(%s, s) ',args(1:end-2));
            s = sprintf('%s %s(%s)./(%s)%s',argF,sint,snum,sden,sdel);
            
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
            else
                ARG = {};
            end
            
            %%% Main
            h = obj.templates(ishow).show(ARG{:});
            
            % plot nominal
            if strcmp(opt,'nom')
                if isempty(obj.nominal)
                    disp('no nominal exists.');
                else
                    %plot(real(obj.nominal),imag(obj.nominal));
                    obj.nominal.show('k');
                end
            elseif ~strcmp(opt,'nonom')
                error('mode options: ''nom'' | ''nonom''')
            end
            
            if nargout==1
                varargout{1}=h;
            end
        end
        function varargout=cases(obj,par,w)
            %CASES returns the template points for given parametric cases
            %It does not plot anything!
            %
            %  Usage:
            %  [T,w] = CASES(obj,par,w)   computes a template T in nichols
            %  foramt for given parameters par, and frequencies w
            %
            %  Inputs:
            %  par      array with each column a different parameter case;
            %           default is the eniter grid
            %  w        frequency vector; default is the frequency vector
            %           of nominal case (if avialble) or logspace(-2,2,50)
            
            if nargin<3, w=[]; end
            if nargin<2, par=[]; end
            
            if isempty(w) && isempty(obj.nominal)
                w = logspace(-2,2,50);
            elseif isempty(w)
                w = obj.nominal.frequency;
            end
            
            w = reshape(w,[],1); % make sure w is a column vector.
            
            f = obj.qplant2func();
            if isempty(par)
                pgrid = grid(obj.pars,[],0);
            elseif size(par,1) ~= length(obj.pars)
                error('par must have %i raws', length(obj.pars));
            else
                pgrid = par;
            end
            p = qplant.pack(pgrid);
            
            varargout{1} = obj.funcval(f,p,w*1j); % response in Nichols format
            if nargout==2
                varargout{2} = w;
            elseif nargout>2
                error('too many outputs!')
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
            if isempty(obj.templates), error('No templates are avilable'); end
            
            Freq = [obj.templates.frequency];
            T = obj.templates(Freq==w);
            if isempty(T), error('requested frequency is not in template'); end
            if isempty(idx)
                tpl = T.template;
                par = T.parameters;
            else
                [tpl,par] = get(T,idx);
            end
            
        end
        function sys = mtimes(A,B)
           %MTIMES series connection of a qplant object
           %    Same as QPLANT/SERIES
           sys = qsys({A,B},'B{1}*B{2}');
        end
        function sys = series(A,B)
            %SERIES series connection of a qplant object
            sys = qsys({A,B},'B{1}*B{2}');
        end
        function sys = feedback(A,B)
            %FEEDBACK connect qplant objects by negative feedback
            sys = qsys({A,B},'1/(1+B{1}*B{2})');
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
                o = length(p)-1;    % highest order
                for k=1:(o+1)
                    s = sprintf('%s %g*s.^%i + ',s,p(k),o);
                    o=o-1;
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
            %case C and for S. The output is in unwrapped Nichols format.
            
            c{end} = s;
            
            %this part is not needed for Matlab 2016a or later
            for k=1:length(c)-1
                c{k} = repmat(c{k},[size(s,1) 1]);
            end
            c{end} = repmat(s,[1 size(c{1},2)]);
            
            nyq  = f(c{:});
            T=c2n(nyq,'unwrap');
        end
        [T,Qpar] = adgrid(trf,s,qpar,Tacc,n,plot_on)
    end
    
end

