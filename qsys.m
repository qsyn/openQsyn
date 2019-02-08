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
                blk = obj.blocks{k};
                if isa(blk,'qplant')
                    nq = nq+1;
                elseif ~any([isa(blk,'qfr') isa(blk,'lti') (isnumeric(blk) && isreal(blk))])
                    error('qsys accepts only qplants, qfr, lti, and real scalar blocks')
                end
            end
            if nq==0
                error('qsys requires at least 1 block element to be a qplant object')
            end
            obj.nqplant = nq;
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
            %       shophase    show phase plot: 0 | 1 (def)
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
            
            if ~iscell(par), par={par}; end % force cell
            
            N = length(obj.blocks); % number of blocks
            T = cell(obj.nqplant,1);
            kp=0;
            for k=1:N
                blk = obj.blocks{k};
                if isa(blk,'qplant')
                    kp = kp+1;
                    T{k} = blk.cases(par{kp},w);
                elseif isnumeric(blk)
                    T{k} = blk;
                else
                    T{k} = freqresp(blk,w);
                end
            end
            f = qsys2func(obj);
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
    end
end

