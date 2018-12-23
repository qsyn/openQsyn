classdef qftpar
    %QFTPAR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name        % name
        nominal     % Nominal value
        lower       % lower bound
        upper       % upper bound
        grid        % number of grid points 
    end
    
    methods
        function par = qftpar(name,nominal,lower,upper,grid)
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


