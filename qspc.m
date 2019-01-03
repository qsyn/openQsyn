classdef qspc
    %QSPC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        frequency
        upper
        lower
    end
    
    methods
        function obj = qspc(name,w,upper,lower)
            %SPC Construct an instance of this class
            %   Detailed explanation goes here
            if nargin<4, lower=[]; end

            % check inputs
            if ~ischar(name), error('qpar name must be a charecter array'); end
            if ~isnumeric(w), error('frqueincy must be a numeric scalar'); end
            if ~isnumeric(upper), error('upper bound must be a numeric scalar'); end
            if ~isnumeric(lower), error('lower bound must be a numeric scalar'); end
            if any(w<0) || min(size(w)>2)
                error('w must be a vector of positive values'); 
            end
            
            if isscalar(upper) 
                upper = upper*ones(1,length(w));
            elseif  length(upper)~=length(w)
                error('upper spec must be scalar or a vector of same length as w');
            end
            
            if ~isempty(lower)
                if isscalar(lower)
                    lower = lower*ones(1,length(w));
                elseif length(lower)~=length(w)
                    error('lower spec must be scalar or a vector of same length as w');
                end
            end

            obj.name = name;
            obj.frequency = w;
            obj.upper = upper;
            obj.lower = lower;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

