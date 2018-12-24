classdef qpoly < matlab.mixin.CustomDisplay
    %QPOLY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        a
        n
        pars
    end
    
    methods

        function obj = qpoly(varargin)
            %QPOLY Construct an instance of this class
            %   Detailed explanation goes here
            if nargin==0, error('qpoly requires at least one argument'); end
            obj.a = cell(nargin);
            Pars = [];
            for k=1:nargin
                switch class(varargin{k})
                    case 'qpar'
                        obj.a{k} = varargin{k};
                        Pars = vertcat(Pars ,varargin{k});
                    case 'qexpression'
                        obj.a{k} = varargin{k};
                        Pars = vertcat(Pars,varargin{k}.pars);
                    case 'double'
                        obj.a{k} = varargin{k};
                    otherwise
                         error('cannot concatenate qpar and whatever you have there')
                end 
            end
            obj.pars = unique(Pars);
            obj.n = length(obj.a)-1;
        end
        function val = nom(obj)
            val = zeros(1,length(obj.a));
            for k=1:length(obj.a)
                if isnumeric(obj.a{k})
                    val(k) = obj.a{k};
                else 
                    val(k) = obj.a{k}.nom;
                end
            end     
        end
        function h = qpoly2func(obj)
            s = '';
            for k=1:obj.n+1
                if isa(obj.a{k},'qexpression')
                    hk = qexp2func(obj.a{k});
                    sk = func2str(hk);
                    idx = strfind(sk,')');  
                    sk = sk(idx+1:end);
                elseif  isa(obj.a{k},'qpar')
                    sk = obj.a{k}.name;
                else
                    sk = num2str(obj.a{k});
                end
                s = strcat(s,sprintf(' (%s).*s.^%i + ',sk,obj.n+1-k));
                
            end
            args = sprintf('%s, ',obj.pars.name);
            argF = sprintf('@(%s,s) ',args(1:end-2));
            s = s(2:end-2);
            h = str2func([argF s]);
        end

    end
    
    methods (Access = protected)
        
        function header = getHeader(obj)
            if ~isscalar(obj)
                header = getHeader@matlab.mixin.CustomDisplay(obj);
            else
                headerStr = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
                headerStr = [headerStr,' with coefficients'];
                header = sprintf('%s\n',headerStr);
            end
        end
        function propgrp = getPropertyGroups(obj)
            if ~isscalar(obj)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);                
            else
                s = cell(length(obj.a),1);
                c = cell(length(obj.a),1);
                for k=obj.n+1:-1:1
                    s{k} = sprintf('s%i',k-1);
                    if isa(obj.a{k},'qpar')
                        c{k} = obj.a{k}.name;
                    elseif isa(obj.a{k},'qexpression')
                        c{k} = obj.a{k}.expression;
                    else
                        c{k} = obj.a{k};
                    end
                end
                propList = cell2struct(c,flipud(s));  

                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        end
        
    end
end

