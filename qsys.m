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
        function bodcases(obj,par,w,opt)
            %BODCASES system frequency domain simulation for user selected
            %cases on Bode plot
            %
            %   BODCASES(sys)   plots bode for all cases given by the 
            %   parameters of the plants in the qsys object sys
            %
            %   BODCASES(QPLANT,W,PAR)   specify frequency and cases set
            %
            %   BODCASES(QPLANT,W,OPT)   specify wihch part to plot: 
            %                            'mag' | 'phase' | 'magphase' (def)
            %
            %   
            
            % TO DO: something with OPT
            
            if nargin<4, opt=[]; end                    
            if nargin<3, w=[]; end
            if nargin<2, par=[]; end    
            if isempty(opt), opt = 'magphase'; end    
            
            [res,w] = cases(obj,par,w);
            col = lines(size(res,2));
            qtpl.bodeplotter(res.',w,opt,col);           
                
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

