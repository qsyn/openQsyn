classdef qsys
    %QSYS is a class used to decribe a system composed of at least one qplant
    %It may include additional blocks such as LTI objects and additional qplant
    %objects
    
    properties
        blocks          % cell array containing all blocks
        connections     % describe connections between blocks
        nqplant         % number of qplant blocks
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
                elseif isa(blk,'qsys')
                    nq = nq + blk.nqplant;
                elseif ~any([isa(blk,'qfr') isa(blk,'lti') isa(blk,'qctrl') (isnumeric(blk) && isreal(blk))])
                    error('qsys accepts only qplants, qctrl, qfr, lti, and real scalar blocks')
                end
            end
            if nq==0
                error('qsys requires at least 1 block element to be a qplant object')
            end
            obj.nqplant = nq;
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
        function varargout = step(varargin)
            % STEP is used to calculate and plot the step response of a
            % QSYS object for a given set of parameters and a time vector.
            %  Usage: 
            %           y = STEP(obj,pars,t)   computes the m step
            %           responses of obj for the m different cases of the n
            %           uncertain parameters given by the nXm matrix
            %           "pars", for the time vector t. 
            %
            %           STEP(obj,pars,t)        plots the responses
            %
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
            %
            %           t       a time vector for the simulation, i.e.
            %                   t=0:0.1:30. If no time vector is supplied,
            %                   the default is 0:0.1:10
            defaultT=0:0.1:10;
            nomPars=1;
            switch nargin
                case 1
                    obj=varargin{1};
                    assert(isa(obj,'qsys'),'Must input a qsys object!')
                    t=defaultT;
                    pars=nomPars;
                case 2
                    obj=varargin{1};
                    assert(isa(obj,'qsys'),'Must input a qsys object!')
                    if (isrow(varargin{2}))
                        t=varargin{2};
                        pars=nomPars;
                    else
                        pars=varargin{2};
                        t=defaultT;
                    end
                case 3
                    obj=varargin{1};
                    pars=varargin{2};
                    t=varargin{3};
            end
            [~,nCases]=size(pars);
            y=zeros(length(t),nCases);
            for ii=1:nCases
                [sysNum,sysDen]=CoeffExtract(obj,pars(:,ii));
                [a,b,c,d] = qplant.Local_tf2ss(sysNum,sysDen);
                x0=zeros(length(a),1);
                [t x] = ode45(@(t,x) odeFun(t,x,a,b), t, x0);
                y(:,ii) = c*x'+d;
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
        function [sysNum,sysDen]=CoeffExtract(obj,pars)
            % COEFFEXTRACT is used to extract the numerator/denominator
            % coefficient vectors of the form
            %               NUM(s)
            %       H(s) = -------
            %               DEN(s)
            % from a QSYS object for time domain simulations.
            %  Usage: 
            %           [SYSNUM,SYSDEN] = COEFFEXTRACT(obj,pars)
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
            connection=obj.connections;
            if isa(obj.blocks{1},'qsys')
                [numBlk1,denBlk1]=CoeffExtract(obj.blocks{1},pars);
            end
            if isa(obj.blocks{2},'qsys')
                [numBlk2,denBlk2]=CoeffExtract(obj.blocks{2},pars);
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
            if (Bool2&&Bool1)
                indInt=find(denBlk2,1,'last');
                indDif=find(numBlk2,1,'last');
                gain=numBlk2(indDif)/denBlk2(indInt);

                Blk2=qctrl(roots(numBlk2),roots(denBlk2),gain);
                Blk1=qplant(numBlk1,denBlk1);
                obj=qsys({Blk1,Blk2},connection);
            end

            [sysNum,sysDen]=ReducedSolver(obj,pars);
        end
        function [sysNum,sysDen]=ReducedSolver(obj,pars)
            % REDUCEDSOLVER solves the base-case of a QSYS object involving
            % a single QPLANT object and another QCTRL/TF/ZPK/DOUBLE
            % object. It outputs the numerator/denominator data of the form
            %               NUM(s)
            %       H(s) = -------
            %               DEN(s)
            % for the interconnected object.            
            %  Usage: 
            %           [SYSNUM,SYSDEN] = REDUCEDSOLVER(obj,pars)
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
            assert((isa(obj,'qsys') && isa(pars,'numeric')),'Input must be (qsys,numeric)')

            feedCon='1/(1+B{1}*B{2})';
            serCon='B{1}*B{2}';
            for k=1:length(obj.blocks)
                blk = obj.blocks{k};
                switch class(blk)
                    case 'zpk'
                        num=blk.K*poly(blk.Z{:}); %maybe need to fix orders
                        den=poly(blk.P{:});
                    case 'qctrl'
                          [num,den] = tfdata(blk);
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
            numProd=conv(numP,num);
            denProd=conv(denP,den);
            padSize=length(denProd)-length(numProd);
            numProd=[zeros(1,padSize),numProd];

            switch obj.connections
                case feedCon
                    sysNum=denProd;
                    sysDen=denProd+numProd;
                case serCon
                    sysNum=numProd;
                    sysDen=denProd;
            end
        end
        
    end
end

