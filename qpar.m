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
        function obj = uplus(A)
            obj = qexpression([],A,'+');
        end
        function obj = minus(A,B)
            if isnumeric(B) || isa(B,'qpar') ||  isa(B,'qexpression') 
                obj = qexpression(A,B,'-');
            else
                error('undefined method');
            end
        end
        function obj = uminus(A)
            obj = qexpression([],A,'-');
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
        function p = grid(obj,rnd,cases)
            %GRID computes a gird of N vectors
            % 
            %   p = grid(obj,rnd,cases)
            %   
            %   Inputs 
            % 
            %   rnd     boolian scalar to specify if grid is random. def
            %           option is 0 (not random). 
            %
            %   cases   intiger specifing the number of cases 
            %
            %%% Adapted from NDGRID
            [N,M] = size(obj);
            if M~=1, error('par.gird only accepts a column vector'); end            
            if nargin<3, cases = [obj.cases]; end
            if nargin<2, rnd=0; end
            if length(cases)==1 && N>1, cases=repmat(cases,1,N); end
            p = zeros(N,prod(cases));
            for k=1:N
                if rnd==0
                    a = linspace(obj(k),cases(k));
                else %rnd==1        
                    a = obj(k).lower + (obj(k).lower+obj(k).upper)*rand(1,cases(k));
                end
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
        function p = sample(obj,N)
            %RNDGRID generates N random samples of the parameter cases
            %   generates a grid 
            [n,m] = size(obj);
            if m~=1, error('par.rndgird only accepts a column vector'); end
            
            p = zeros(n,N);
            
            for k=1:n
                p(k,:) = obj(k).lower + (obj(k).lower+obj(k).upper)*rand(1,N);
            end
                        
        end
        function I = ismember(A,B)
            %ISMEMBER returns a vector of logical indices positive for
            %every element of parameter set A that is a member of parameter
            %set B.
            Anames = {A(:).name};
            Bnames = {B(:).name};
            I = ismember(Anames,Bnames);
        end
    end

end


