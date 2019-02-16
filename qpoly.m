classdef qpoly < matlab.mixin.CustomDisplay
%QPOLY a class for polynomials containing qpar and/or qexpression objects 
%   
%   A QPOLY object mayt be created by one of the following methods:
%   
%   p = [a b c ...]   constructs a QPOLY object with coefficienct a,b,c,...  
%   Each can be a qpar, qexpression, numeric scalar, 
%
%   p = QPOLY(a,b,c,...)   equivalent use
%
%   p = QPOLY(SYMS,PARS)   constructs a QPOLY object from a symbolic 
%   expression string SYMS and a qpar row vector PARS 
%
%For a list QPOLY methods type: methods qpoly
%For help on a specific methods type: help qpoly/<method>
%
%See also: QPAR, QEXPRESSION, QPLANT

    properties
        a
        n
        pars
    end
    
    methods
        function obj = qpoly(varargin)
            %QPOLY Construct an instance of the QPOLY class
            %   
            % obj=QPOLY(a,b,c,...)   constructs a QPOLY object with
            % coefficienct given by a,b,c,...
            %
            % obj=QPOLY(SYMS,PARS)   constructs a QPOLY object from a
            % symbolic expression SYMS containning the parameters in PARS
           
            if nargin<2, error('qpoly requires at least two argument'); end
            
            if isa(varargin{1},'sym') && isa(varargin{2},'qpar')
                obj = qpoly.sym2qpoly(varargin{1},varargin{2});
                return
            elseif isa(varargin{1},'sym') && ~isa(varargin{2},'qpar')
                error('if input 1 is a symbolic expression, input 2 must be a QPAR array');
            end
                
            obj.a = cell(1,nargin);
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
            %NOM computes the nominal polynomial
            val = zeros(1,length(obj.a));
            for k=1:length(obj.a)
                if isnumeric(obj.a{k})
                    val(k) = obj.a{k};
                else 
                    val(k) = obj.a{k}.nom;
                end
            end     
        end
        function val = cases(obj,par)
            %POLYCASE returns the polynom coeffiecients for a given
            %parameter case
            if nargin<2, par=[]; end
            
            val = zeros(1,length(obj.a));
            for k=1:length(obj.a)
                if isnumeric(obj.a{k})
                    val(k) = obj.a{k};
                elseif isempty(par)
                    val(k) = obj.a{k}.nom;
                else
                    ip = ismember(obj.pars,obj.a{k}.pars);
                    par1 = par(ip);
                    val(k) = obj.a{k}.cases(par1);
                end
            end     
        end
        function h = qpoly2func(obj)
            %QPOLY2FUNC return an handle to a function object f@(p1,p2,...pn,s) 
            %with p1,...,pn corresponding to uncertain parameters.
            
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
        function varargout = coeffs(obj)
           %COEFFS returns the qpoly coefficients
           %
           %    [an,...,a1,a0] = coeffs(p)  returns all coefficeints by order
           % 
           for k=1:nargout
               varargout{k} = obj.a{k};
           end
        end
    end
    
    methods(Access = protected)
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
    
    methods(Static)
        function p = sym2qpoly(S,pars)
            %SYM2POLY creates qpoly from symbolic expression and qpar array
            %
            %Inputs:
            %1) S       sybmolic expression generated by theSymbolic toolbox
            %2) pars    column vector of qpar objects
            
            c=coeffs(S,'s');
            e = cell(1,length(c));
            for k=1:length(c)
                e{end-k+1} = qexpression(char(c(k)),pars);
            end
            p = [e{:}];
        end
    end
end

