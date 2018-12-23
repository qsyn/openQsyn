classdef qpar
    %QPAR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name        % name
        nominal     % Nominal value
        lower       % lower bound
        upper       % upper bound
        grid        % number of grid points 
    end
    
    methods
        function par = qpar(name,nominal,lower,upper,grid)
            %QPAR Construct an instance of this class
            %   Detailed explanation goes here
            
            % check input types
            if nargin~=4, error('qpar constructor demands 4 inputs'); end
            if ~ischar(name), error('qpar name must be a charecter array'); end
            if ~isnumeric(ubnd), error('upper bound must be a numeric scalar'); end
            if ~isnumeric(lbnd), error('lower bound must be a numeric scalar'); end
            if ~isnumeric(nom), error('nominal value must be a numeric scalar'); end
            if any(size(ubnd)~=1), error('upper bound must be a numeric scalar'); end
            if any(size(lbnd)~=1), error('lower bound must be a numeric scalar'); end
            if any(size(nom)~=1), error('nominal value must be a numeric scalar'); end
            
            % check bounds
            if (nom > ubnd) || (nom <lbnd), error('nominal value must be between bounds'); end
            
            % assign properties    
            par.name = char(name);
            par.nominal = double(nominal);
            par.lower = double(lower);
            par.upper = double(upper);
            if nargin==5
                par.grid = int32(grid);
            else
                par.grid = 3; 
                disp('3 grid points are selected as default');
            end
        end
        function obj = plus(A,B)
            if isnumeric(B) || isa(B,'qftpar') ||  isa(B,'qftexpression')
                obj = qftexpression(A,B,'+');
            else
                error('undefined method');
            end
        end
        function obj = minus(A,B)
            if isnumeric(B) || isa(B,'qftpar') ||  isa(B,'qftexpression') 
                obj = qftexpression(A,B,'-');
            else
                error('undefined method');
            end
        end
        function obj = mtimes(A,B)
            if isnumeric(B) || isa(B,'qftpar') ||  isa(B,'qftexpression')
                obj = qftexpression(A,B,'*');
            else
                error('undefined method');
            end
        end
        function obj = mrdivide(A,B)
            if isnumeric(B) || isa(B,'qftpar') ||  isa(B,'qftexpression')
                obj = qftexpression(A,B,'/');
            else
                error('undefined method');
            end
        end
        function par = unique(par)
            names = {par.name};
            [~,ia] = unique(names);
            par = par(ia);
        end      
    end
    
end


