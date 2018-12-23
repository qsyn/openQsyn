classdef qpar
    %QPAR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name
        LowerBnd
        UpperBnd
        Nominal
    end
    
    methods
        function obj = qpar(name,lbnd,ubnd,nom)
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
            obj.Name = char(name);
            obj.LowerBnd = double(lbnd);
            obj.UpperBnd = double(ubnd);
            obj.Nominal = double(nom);
        end
        
        function uobj = unique(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Names = {[]};
            pars
            for k=1:length(obj)
                
                Names{k}
                
            end
            
        end
    end
end

