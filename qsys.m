classdef qsys
    %QSYS is a class used to decribe a system composed of at least one qplant
    %It may include additional blocks such as LTI objects and additional qplant
    %objects
    
    properties
        blocks          % cell array containing all blocks
        connections     % describe connections between blocks
        nqplant         % number of qplant blocks
        delay           % maximal internal delay
    end
    
    methods
        function obj = qsys(blocks,connections)
            %QSYS Construct an instance of this class
            %   Detailed explanation goes here
            p = inputParser;
            valB = @(x) validateattributes(x,{'cell','qplant'},{'nonempty'});
            % || validateattributes(x,{'qplant'},{'scalar'});
            valC = @(x) validateattributes(x,{'char'},{'nonempty'});
            addRequired(p,'blocks',valB);
            addRequired(p,'connections',valC);
            parse(p,blocks,connections)
            if iscell(p.Results.blocks)
                obj.blocks = p.Results.blocks;
            else
                obj.blocks = {p.Results.blocks}; % force into cell
            end
            obj.connections = p.Results.connections;
            
            nq=0; %qplant counter
            for k=1:length(obj.blocks)
                blk = obj.blocks{k}; %correct the index
                if isa(blk,'qplant')
                    nq = nq+1;
                    if ~isempty(blk.delay)
                        if (isempty(obj.delay) || blk.delay>obj.delay)
                            obj.delay=blk.delay;
                        end
                    end
                elseif isa(blk,'qsys')
                    nq = nq + blk.nqplant;
                    if ~isempty(blk.delay)
                        if (isempty(obj.delay) || blk.delay>obj.delay)
                            obj.delay=blk.delay;
                        end
                    end
                elseif ~any([isa(blk,'qfr') isa(blk,'lti') isa(blk,'qctrl') (isnumeric(blk) && isreal(blk))])
                    error('qsys accepts only qplants, qctrl, qfr, lti, and real scalar blocks')
                end
            end
            if nq==0
                error('qsys requires at least 1 block element to be a qplant object')
            end
            obj.nqplant = nq;
            if isempty(obj.delay)
                obj.delay=0;
            end
        end
        function obj = series(A,B)
            %SERIES series connection of a qsys object
            obj = qsys({A,B},'B{1}*B{2}');
        end
        function bodcases(obj,varargin)
            %BODCASES system frequency response for selected cases on Bode plot
            %
            %   BODCASES(sys)   plots Bode for all cases given by the
            %   parameters in all plants found within the system
            %
            %   BODCASES(sys,par,w)   specify frequency and cases set
            %
            %   BODCASES(sys,par,w,parameter,value)   specify additional options
            %   using parameter/value pairs
            %
            %   qplant.BODCASES(...)    alternative usage
            %
            %
            %   Inputs (Optional):
            %       par     array with each column a different parameter case;
            %               default is the eniter grid
            %       w       frequency vector. default is the frequency vector of
            %               nominal case, if available, or logspace(-2,2,50)
            %
            %   Additional Options:
            %       color       specifiy color: RGB array
            %       shownom     show nominal in bold dashed line: 0 (def) | 1
            %       showmag     show magnitude plot: 0 | 1 (def)
            %       showphase    show phase plot: 0 | 1 (def)
            %
            %   See also: qplant/bodcases
            
            % parse inputs
            p = inputParser;
            addOptional(p,'par',[],@(x) validateattributes(x,{'numeric'},{'2d'},1))
            addOptional(p,'w',[],@(x) validateattributes(x,{'numeric'},{'positive'},1))
            addParameter(p,'color',[],@(x) validateattributes(x,{'numeric'},{'ncols', 3}))
            addParameter(p,'showmag',1,@(x) validateattributes(x,{'numeric'},{'scalar','binary'}))
            addParameter(p,'showphase',1,@(x) validateattributes(x,{'numeric'},{'scalar','binary'}))
            addParameter(p,'shownom',0,@(x) validateattributes(x,{'numeric'},{'scalar','binary'}))
            parse(p,varargin{:})
            
            par = p.Results.par;
            w = p.Results.w;
            col = p.Results.color;
            shownom = p.Results.shownom;
            
            if p.Results.showmag && p.Results.showphase
                opt = 'magphase';
            elseif p.Results.showmag
                opt = 'mag';
            else
                opt = 'phase';
            end
            
            [res,w] = cases(obj,par,w);     % compute cases
            N = size(res,2);
            if isempty(col)                 % pick colors
                col = lines(size(res,2));
            elseif size(col,1) < N
                col = repmat(col,ceil(N/size(col,1)),1);
            end
            
            linespec = struct('width',1,'style','-');
            bodeplotter(res.',w,opt,col,linespec);   % plot the Bode for all cases
            
            if  shownom
                linespec.width = 3;
                linespec.style = '--';
                nompar = obj.pars.nom;
                nomres = cases(obj,nompar,w);
                bodeplotter(nomres.',w,opt,[0 0 0],linespec); % plot the Bode for nominal case
            end
            
        end
        function niccases(obj,varargin)
            %NICCASES system frequency response for selected cases on Nichols chart
            %
            %   NICCASES(sys)   plots Nichols chart for all cases given by the
            %   plant parameters
            %
            %   NICCASES(sys,par,w)   specify frequency and cases set
            %
            %   NICCASES(sys,par,w,parameter,value)   specify additional options
            %   using parameter/value pairs
            %
            %   qplant.NICCASES(...)    alternative usage
            %
            %
            %   Inputs (Optional):
            %       par     array with each column a different parameter case;
            %               default is the eniter grid
            %       w       frequency vector. default is the frequency vector of
            %               nominal case, if available, or logspace(-2,2,50)
            %
            %   Additional Options:
            %       color       RGB array
            %       shownom     show nominal in bold dashed line 0 (def) | 1
            %
            %   See also: qplant/niccases
            
            p = inputParser;
            addOptional(p,'par',[],@(x) validateattributes(x,{'numeric'},{'2d'},1))
            addOptional(p,'w',[],@(x) validateattributes(x,{'numeric'},{'positive'},1))
            addParameter(p,'color',[],@(x) validateattributes(x,{'numeric'},{'ncols', 3}))
            addParameter(p,'shownom',0,@(x) validateattributes(x,{'numeric'},{'scalar','binary'}))
            parse(p,varargin{:})
            
            par = p.Results.par;
            w = p.Results.w;
            col = p.Results.color;
            shownom = p.Results.shownom;
            
            [res,w] = cases(obj,par,w);     % compute cases
            N = size(res,2);
            if isempty(col)                 % pick colors
                col = lines(size(res,2));
            elseif size(col,1) < N
                col = repmat(col,ceil(N/size(col,1)),1);
            end
                       
            linespec = struct('width',1,'style','-');
            nicholsplotter(res.',-w,col,linespec);  % plot the Bode for all cases.
                                                    % the (-w) is to get the frequency
                                                    % data tip dispalyed right.
            
            if  shownom
                hold on
                linespec.width = 3;
                linespec.style = '--';
                nompar = obj.pars.nom;
                nomres = cases(obj,nompar,w);
                nicholsplotter(nomres.',w,[0 0 0],linespec);   % plot the Nichols for nominal case
            end
            
        end
        function varargout = cases(obj,par,w)
            %CASES returns the template points for given parametric cases
            %It does not plot anything!
            %
            %  Usage:
            %  [T,w] = CASES(obj,par,w)   computes a template T in nichols
            %  foramt for given parameters par, and frequencies w
            %
            %  Inputs:
            %  par      array with each column a different parameter case;
            %           default is the eniter grid. for multiple qplant
            %           inset pars for each plant in cells: par={p1,p2,...}
            %  w        vector of frequencies, def = logspace(-2,3,200)
            
            if nargin<3, w=[]; end
            if nargin<2, par=[]; end
            
            if isempty(w), w = logspace(-2,3,200); end
            w = reshape(w,[],1); % make sure w is a column vector.
            
            if isempty(par), par=cell(obj.nqplant,1);  end  % no par specified
            if ~iscell(par), par={par}; end                 % force cell
            
            N = length(obj.blocks); % number of blocks
            T = cell(obj.nqplant,1);
            kp=0;
            for k=1:N
                blk = obj.blocks{k};
                if (isa(blk,'qplant') || isa(blk,'qsys'))
                    kp = kp+1;
                    if isscalar(par)
                        T{k} = n2c(blk.cases(par{1},w));
                    else
                        T{k} = n2c(blk.cases(par{kp},w));
                    end
                else
                    T{k} = n2c(nicresp(blk,w)).';
                end
            end
            f = qsys2func(obj);
            %[nrows,ncols] = cellfun(@size,T); % continue here...
            
            %for k=1:N
                
            %end
            varargout{1} = c2n(f(T),-180);
            if nargout==2
                varargout{2} = w;
            elseif nargout>2
                error('too many outputs!')
            end
            
        end
        function f = qsys2func(obj)
            %QSYS2FUNC creates function handle according to coonection
            argF = '@(B) ';
            exp = replace(obj.connections,{'*','/','^'},{'.*','./','.^'});
            f = str2func([argF exp]);
        end
        function sys = feedback(A,B)
            %FEEDBACK connects qplant or qsys objects by negative feedback
            sys = qsys({A,B},'1/(1+B{1}*B{2})');
        end
        function varargout = step(obj,varargin)
        % STEP is used to calculate and plot the step response of a
        % QSYS object for a given set of parameters and a time vector.
        %  Usage: 
        %           y = STEP(obj,'Pars',pars,'Time',t)   computes the m step
        %           responses of obj for the m different cases of the n
        %           uncertain parameters given by the nXm matrix
        %           "pars", for the time vector t. If a delay is present
        %           everything is discretized and the delay is absorbed.
        %
        %           STEP(obj,'Pars',pars,'Time',t)       plots the responses
        %
        %  Inputs:
        %           obj       a QSYS object representing an
        %                     interconnection of systems.
        %
        %           varargin  name-parameter pairs to be parsed by the function.
        %
        %           'Pars'      an array with n rows, each corresponds to
        %                     a QPAR object, and m colums each
        %                     corresponds to a test case to be simulated.
        %                     For example if the uncertain plant inside
        %                     is denoted P with P.pars.name a,k,wn,z,
        %                     then pars would have 4 rows. If pars is
        %                     empty, the nominal values will be
        %                     simulated.
        %
        %           'Time'    a time vector for the simulation, i.e.
        %                     t=0:0.1:30. If no time vector is supplied,
        %                     the default is 0:0.1:10
        %
        %           'Ts'      a scalar representing desired sample time for
        %                     discretization. Currently deals with delays, future
        %                     support to computer-controlled controllers.
                % Input parsing
                    defaultT=0:0.1:10;
                    nomPars=1;
                    defaultTs=0;
                    h=0;
                    nu=20;
                    p = inputParser;
                    validObject = @(x) isa(x,'qsys');
                    validTs= @(x) isnumeric(x) && isscalar(x);
                    addRequired(p,'obj',validObject);
                    parse(p,obj);
                    if obj.delay
                        h=obj.delay;
                        defaultTs=h/nu;
                        validTs= @(x) isnumeric(x) && isscalar(x) && ~mod(h,x);
                    end
                    addParameter(p,'Pars',nomPars,@ismatrix);
                    addParameter(p,'Time',defaultT,@isvector);
                    addParameter(p,'Ts',defaultTs,validTs);
                    parse(p,obj,varargin{:});

                    Ts=p.Results.Ts;
                    pars=p.Results.Pars;
                    t=p.Results.Time;
                    obj=p.Results.obj;
                    [~,nCases]=size(pars);
                    % If there is delay discretize, else don't
                    if h == 0
                        y=zeros(length(t),nCases);
                        for ii=1:nCases
                            [sysNum,sysDen]=CoeffExtract(obj,pars(:,ii),h,Ts);
                            [a,b,c,d] = Local_tf2ss(sysNum,sysDen);
                            x0=zeros(length(a),1);
                            [t x] = ode45(@(t,x) odeFun(t,x,a,b), t, x0);
                            y(:,ii) = c*x'+d;
                        end
                    else
                        t=0:Ts:max(t);
                        for ii=1:nCases
                            [sysNum,sysDen]=CoeffExtract(obj,pars(:,ii),h,Ts);
                            [a,b,c,d] = Local_tf2ss(sysNum,sysDen);
                            clear x
                            x0=zeros(length(a),1);
                            x(:,1) = a*x0 + b.*1;
                            y(ii,1) = c*x(:,1);
                        % x(k + 1) = A*x(k) + B*u(k)
                            for jj = 2:length(t)
                                x(:,jj) = a*x(:,jj-1) + b*1;
                                y(ii,jj)=c*x(:,jj)+d*1;
                            end
                        end
                    end
                    col=lines(nCases);
                    linespec = struct('width',1,'style','-');
                    figure
                        hold on
                        set(gca, 'ColorOrder', col, 'NextPlot', 'add')
                        plot(t,y,'linewidth',linespec.width,'linestyle',linespec.style);
                        xlabel('Time [s]');
                        ylabel('Amplitude');
                        title('Step Response')
                        axis tight
                        hold off
                        if nargout
                            varargout{1}=y;
                        end
                    function dxdt = odeFun(t,x,A,B)
                    u = 1; % for step input
                    dxdt = A*x+B; % simply write the equation
                    end 
        end
        function varargout = lsim(obj,varargin)
        % LSIM is used to calculate and plot the linear response of a
        % QSYS object for a given set of parameters, an input and a time vector.
        %  Usage: 
        %           y = LSIM(obj,'Pars',pars,'Time',t,'u',u)   computes the m linear responses of obj
        %           to a forcing signal u. There are m different cases of the n
        %           uncertain parameters given by the nXm matrix "pars", for the
        %           time vector t. If a delay is present everything is discretized
        %           and the delay is absorbed into the parameters. If nargout = 0,
        %           simply plot the responses.
        %
        %  Inputs:
        %           obj       a QSYS object representing an
        %                     interconnection of systems.
        %           
        %           varargin  name-parameter pairs to be parsed by the function.
        %
        %           'Pars'    an array with n rows, each corresponds to
        %                     a QPAR object, and m colums each
        %                     corresponds to a test case to be simulated.
        %                     For example if the uncertain plant inside
        %                     is denoted P with P.pars.name a,k,wn,z,
        %                     then pars would have 4 rows. If pars is
        %                     empty, the nominal values will be
        %                     simulated.
        %
        %           'Time'    a time vector for the simulation, i.e.
        %                     t=0:0.1:30. If no time vector is supplied,
        %                     the default is 0:0.1:10
        %
        %           'u'       a vector representation of a user generated signal.
        %              
        %           'Ts'      a scalar representing desired sample time for
        %                     discretization. Currently deals with delays, future
        %                     support to computer-controlled controllers.
                   %Parse input
                    defaultT=0:0.1:10;
                    nomPars=1;
                    defaultTs=0;
                    h=0;
                    nu=20; % Default integer multiple of delay
                    defaultU=ones(1,length(defaultT));
                    p = inputParser;
                    validObject = @(x) isa(x,'qsys');
                    validTs= @(x) isnumeric(x) && isscalar(x);
                    validU = @(x) isnumeric(x) && isvector(x);
                    addRequired(p,'obj',validObject);
                    parse(p,obj);
                    if obj.delay
                        h=obj.delay;
                        defaultTs=h/nu; 
                        % Default sample time can be adjusted here
                        validTs = @(x) isnumeric(x) && isscalar(x) && ~mod(h,x);
                        %Make sure Ts isinteger multiple of delay
                    end
                    addParameter(p,'Pars',nomPars,@ismatrix);
                    addParameter(p,'Time',defaultT,@isvector);
                    addOptional(p,'Ts',defaultTs,validTs);
                    addParameter(p,'u',defaultU,validU);
                    parse(p,obj,varargin{:});
                    Ts=p.Results.Ts;
                    pars=p.Results.Pars;
                    t=p.Results.Time;
                    obj=p.Results.obj;
                    validateattributes(p.Results.u,{'numeric'},{'size',size(t)});
                    u=p.Results.u;

                % Simulation
                [~,nCases]=size(pars);
                if h==0
                    y=zeros(length(t),nCases);
                    tt=t;
                    for ii=1:nCases
                        [sysNum,sysDen]=CoeffExtract(obj,pars(:,ii),h,Ts);
                        [a,b,c,d] = Local_tf2ss(sysNum,sysDen);
                        x0=zeros(length(a),1);
                        [t x] = ode45(@(t,x) odeFun(t,x,a,b,tt,u), t, x0);
                        y(:,ii) = c*x'+d*u;
                    end
                else
                    t2=0:Ts:max(t);
                    for ii=1:nCases
                        [sysNum,sysDen]=CoeffExtract(obj,pars(:,ii),h,Ts);
                        [a,b,c,d] = Local_tf2ss(sysNum,sysDen);
                        clear x;
                        x0=zeros(length(a),1);
                        x(:,1) = a*x0 + b.*interp1(t,u,t2(1));
                        y(ii,1) = c*x(:,1)+d.*interp1(t,u,t2(1));
                    % x(k + 1) = A*x(k) + B*u(k)
                        for jj = 2:length(t2)
                            x(:,jj) = a*x(:,jj-1) + b.*interp1(t,u,t2(jj));
                            y(ii,jj)=c*x(:,jj)+d.*interp1(t,u,t2(jj));
                        end
                    end
                end
                % Plot
                col=lines(nCases);
                linespec = struct('width',1,'style','-');
                figure
                    hold on
                    set(gca, 'ColorOrder', col, 'NextPlot', 'add')
                    if length(y)>length(t)
                        plot(t2,y,'linewidth',linespec.width,'linestyle',linespec.style);
                    else
                        plot(t,y,'linewidth',linespec.width,'linestyle',linespec.style);
                    end
                    xlabel('Time [s]');
                    ylabel('Amplitude');
                    title('Linear Simulation Result')
                    axis tight
                    hold off
                if nargout
                    varargout{1}=y;
                end
            function dxdt = odeFun(t,x,A,B,tt,u)
                % Interpulate between samples
                intU = interp1(tt,u,t);
                dxdt = A*x+B*intU; % simply write the equation
                end 
            end
        function [sysNum,sysDen]=CoeffExtract(obj,pars,h,Ts)
        % COEFFEXTRACT is used to extract the numerator/denominator
        % coefficient vectors of the form
        %               NUM(s)
        %       H(s) = -------
        %               DEN(s)
        % from a QSYS object for time domain simulations.
        %  Usage: 
        %           [SYSNUM,SYSDEN] = COEFFEXTRACT(obj,pars,h)
        %           recursively finds the basic QSYS interconnection
        %           involving one QPLANT and antoher system and calls
        %           for REDUCEDSOLVER which solves the base case.
        %
        %  Inputs:
        %           obj     a qsys object representing an
        %                   interconnection of systems.
        %
        %           pars    an array with n rows, each corresponds to
        %                   a QPAR object, and m colums each
        %                   corresponds to a test case to be simulated.
        %                   For example if the uncertain plant inside
        %                   is denoted P with P.pars.name a,k,wn,z,
        %                   then pars would have 4 rows. If pars is
        %                   empty, the nominal values will be
        %                   simulated.
        %           h       the maximal delay of the compound qsys object.
        %           Ts      the sample time
            connection=obj.connections;
            %h=obj.delay;
            %nu=40; % Integer multipile of delay
            %Ts=h/nu;
            % Need to fix the Ts as input!
            if isa(obj.blocks{1},'qsys')
                [numBlk1,denBlk1]=CoeffExtract(obj.blocks{1},pars,h,Ts);
            end
            if isa(obj.blocks{2},'qsys')
                [numBlk2,denBlk2]=CoeffExtract(obj.blocks{2},pars,h,Ts);
            end

            Bool1=exist('numBlk1','var');
            Bool2=exist('numBlk2','var');
            if (Bool1&&~Bool2)       
                Blk1=qplant(numBlk1,denBlk1); %%% FINISH THIS!
                obj=qsys({Blk1,obj.blocks{2}},connection);
            end
            if (Bool2&&~Bool1)       
                Blk2=qplant(numBlk2,denBlk2); 
                obj=qsys({obj.blocks{1},Blk2},connection);
            end
            if (Bool2&&Bool1) %Maybe the problem is here
                indInt=find(denBlk2,1,'last');
                indDif=find(numBlk2,1,'last');
                gain=numBlk2(indDif)/denBlk2(indInt);

                Blk2=qctrl(roots(numBlk2),roots(denBlk2),gain);
                if h>0
                    Blk2.sampleTime=Ts;
                end
                Blk1=qplant(numBlk1,denBlk1);
                obj=qsys({Blk1,Blk2},connection);
            end
            [sysNum,sysDen]=ReducedSolver(obj,pars,h,Ts);
        end
        function [sysNum,sysDen]=ReducedSolver(obj,pars,h,Ts)
        % REDUCEDSOLVER solves the base-case of a QSYS object involving
        % a single QPLANT object and another QCTRL/TF/ZPK/DOUBLE
        % object. It outputs the numerator/denominator data of the form
        %               NUM(s)
        %       H(s) = -------
        %               DEN(s)
        % for the interconnected object. If a delay is present the data
        % is discretized and the delay is absorbed as poles at the origin.
        % Currently supports only delays which are integer multiples of Ts.
        % 
        %  Usage: 
        %           [SYSNUM,SYSDEN] = REDUCEDSOLVER(obj,pars,Ts)
        %           finds the basic QSYS interconnection
        %           involving one QPLANT and antoher system. It solves
        %           the base case of COEFFEXTRACT
        %  Inputs:
        %           obj     a QSYS object representing an
        %                   interconnection of systems.
        %
        %           pars    an array with n rows, each corresponds to
        %                   a QPAR object, and m colums each
        %                   corresponds to a test case to be simulated.
        %                   For example if the uncertain plant inside
        %                   is denoted P with P.pars.name a,k,wn,z,
        %                   then pars would have 4 rows. If pars is
        %                   empty, the nominal values will be
        %                   simulated. 
        %           Ts      desired sample time for discretization
        %                   if 0, system is continuous.
        %    assert((isa(obj,'qsys') && isa(pars,'numeric')),'Input must be (qsys,numeric)')
                    feedCon='1/(1+B{1}*B{2})';
                    serCon='B{1}*B{2}';
                    sampleTime=0;
                    for k=1:length(obj.blocks)
                        blk = obj.blocks{k};
                        switch class(blk)
                            case 'zpk'
                                num=blk.K*poly(blk.Z{:}); %maybe need to fix orders
                                den=poly(blk.P{:});
                            case 'qctrl'
                                  [num,den] = tfdata(blk);
                                  sampleTime=blk.sampleTime;
                            case 'tf'
                                num=blk.Numerator{:};
                                den=blk.Denominator{:};
                            case 'double'
                                num=blk;
                                den=1;
                            case 'qplant'
                                Plant=blk;
                        end
                    end
                    if ~isempty(Plant.pars)
                        if isscalar(pars)
                           numP=cases(Plant.num);
                           denP=cases(Plant.den);
                        else
                            [numpars,denpars]=Prase_params(Plant,pars);
                            numP=cases(Plant.num,numpars);
                            denP=cases(Plant.den,denpars);
                        end
                    else
                        numP=Plant.num;
                        denP=Plant.den;
                    end
                    % account for discretization/delays
                    if ~(Ts == 0) && ~(h == 0)
                        % First for the plants
                        if ~isempty(Plant.delay)
                            [num_d,den_d]=AbsorbDelay(numP,denP,Plant.delay,Ts);
                            numP=num_d;
                            denP=den_d;
        %                 else 
        %                     [Phi,Gamma,c,d]=Local_c2d(numP,denP,Ts);
        %                     [num_d,den_d]=Local_ss2tf(Phi,Gamma,c,d,'d');
        %                     numP=num_d;
        %                     denP=den_d;
                        end
                        % Second for the controller (also gain fix)
                        if (length(den)>1 && sampleTime==0)
                            [Phi,Gamma,c,d]=Local_c2d(num,den,Ts);
                            [num_d,den_d]=Local_ss2tf(Phi,Gamma,c,d,'d');
                            K=abs(d-c*(Phi\Gamma));
                            num=K*num_d;
                            den=den_d;
                        end
                    end
                    %

                    switch obj.connections
                        case feedCon
                            numProd=conv(numP,num);
                            denProd=conv(denP,den);
                            padSize=length(denProd)-length(numProd);
                            numProd=[zeros(1,padSize),numProd];
                            sysNum=denProd;
                            sysDen=denProd+numProd;
                        case serCon
                            [numP,den]=qMinreal(numP,den);
                            [num,denP]=qMinreal(num,denP);
                            numProd=conv(numP,num);
                            denProd=conv(denP,den);
                            padSize=length(denProd)-length(numProd);
                            numProd=[zeros(1,padSize),numProd];
                            sysNum=numProd;
                            sysDen=denProd;
                    end
                    %cancel poles and zeros
                     [sysNum,sysDen]=qMinreal(sysNum,sysDen);
        end
    end
end

