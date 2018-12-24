classdef qpar
    %QPAR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name        % name
        nominal     % Nominal value
        lower       % lower bound
        upper       % upper bound
        cases       % number of cases
    end
    
    methods
        function par = qpar(name,nom,lbnd,ubnd,cases)
            %QPAR Construct an instance of this class
            %   Detailed explanation goes here

            if nargin==4 
                cases = 3;
                disp('3 cases are selected as default');
            elseif nargin~=5
                error('unexceptable number of argumetns to qpar')
            end
            
            % check input types
            if ~ischar(name), error('qpar name must be a charecter array'); end
            if ~isnumeric(lbnd), error('lower bound must be a numeric scalar'); end
            if ~isnumeric(ubnd), error('upper bound must be a numeric scalar'); end
            if ~isnumeric(nom), error('nominal value must be a numeric scalar'); end
            if ~isnumeric(cases), error('number of cases must be an integer'); end
            if any(size(ubnd)~=1), error('upper bound must be a numeric scalar'); end
            if any(size(lbnd)~=1), error('lower bound must be a numeric scalar'); end
            if any(size(nom)~=1), error('nominal value must be a numeric scalar'); end
            
            % check bounds
            if (nom > ubnd) || (nom <lbnd), error('nominal value must be between bounds'); end
                       
            % assign properties    
            par.name = char(name);
            par.nominal = double(nom);
            par.lower = double(lbnd);
            par.upper = double(ubnd);
            par.cases = int32(cases);
            
        end
        function obj = plus(A,B)
            if isnumeric(B) || isa(B,'qpar') ||  isa(B,'qexpression')
                obj = qexpression(A,B,'+');
            else
                error('undefined method');
            end
        end
        function obj = minus(A,B)
            if isnumeric(B) || isa(B,'qpar') ||  isa(B,'qexpression') 
                obj = qexpression(A,B,'-');
            else
                error('undefined method');
            end
        end
        function obj = mtimes(A,B)
            if isnumeric(B) || isa(B,'qpar') ||  isa(B,'qexpression')
                obj = qexpression(A,B,'*');
            else
                error('undefined method');
            end
        end
        function obj = mrdivide(A,B)
            if isnumeric(B) || isa(B,'qpar') ||  isa(B,'qexpression')
                obj = qexpression(A,B,'/');
            else
                error('undefined method');
            end
        end
        function obj = unique(par)
            names = {par.name};
            [~,ia] = unique(names);
            obj = par(ia);
        end      
        function obj = horzcat(varargin)
            obj = qpoly(varargin{:});
        end
        function val = nom(par)
            val = par.nominal;
        end
        function pspace = linspace(obj,n)
            if nargin<2, n = obj.cases; end
            pspace = linspace(obj.lower,obj.upper,n);           
        end
        function p = grid(obj,cases)
            %GRID computes a gird of N vectors
            %%% this part adapted from NDGRID
            [N,M] = size(obj);
            if M~=1, error('par.gird only accepts a column vector'); end
            if nargin<2, cases = [obj.cases]; end
            
            p = zeros(N,prod(cases));
            
            for k=1:N
                a = linspace(obj(k),cases(k));
                x = full(a);
                s = ones(1,N);
                s(k) = numel(x);
                x = reshape(x,s);
                s = cases;
                s(k) = 1;
                pk = repmat(x,s);
                p(k,:) = reshape(pk,1,[]);
            end
                        
        end
    end

end


