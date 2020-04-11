classdef qrff
    %QRFF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type
        par1
        par2
    end
    
    methods
        function obj = qrff(type,par1,varargin)
            %QRFF Construct an instance of the qrff class
            %   Detailed explanation goes here
            
            validStrings = {'gain','delay','dc','hf','uns','int'};
            validatestring(type,validStrings,'qrff','type',1)
            p = inputParser;
            p.addRequired('par1',@(x) isscalar(x) && (isnumeric(x) || isa(x,'qpar')));
            p.addOptional('par2',[],@(x) isscalar(x) && (isnumeric(x) || isa(x,'qpar')));
            p.parse(par1,varargin{:});
            
            obj.type = type;
            obj.par1 = p.Results.par1;
            obj.par2 = p.Results.par2;
            
        end
        
        function h = qrff2func(obj)
            %QRFF2FUNC return an handle to a function object f@(p1,p2,...pn,s)
            %with p1,...,pn corresponding to uncertain parameters.
            
            
            Pars = {};
            for k = 1:length(obj)
                if isa(obj(k).par1,'qpar')
                    par1s = obj(k).par1.name;
                    Pars = {Pars{:},par1s};
                else
                    par1s = sprintf('%g',obj(k).par1);
                end
                if isa(obj(k).par2,'qpar')
                    par2s = obj(k).par2.name;
                    Pars = {Pars{:},par2s};
                else
                    par2s = sprintf('%g',obj(k).par2);
                end
                switch obj(k).type
                    case 'gain'
                        Sk = par1s;
                    case 'delay'
                        Sk = sprintf('exp(-s*%s)',par1s);
                    case 'dc'
                        if isempty(obj.par2)
                            Sk = sprintf('(1+s/%s)',par1s);
                        else
                            SdW = sprintf('s./%s',par1s);
                            Sk = sprintf('(1 + 2.*%s.*(%s) + (%s).^2)',par2s,SdW,SdW);
                        end
                    case 'hf'
                        if isempty(obj.par2)
                            Sk = sprintf('(s+%s)',par1s);
                        else
                            wn = par1s;
                            zeta = par2s;
                            Sk = sprintf('(s.^2 + 2.*%s.*%s + %s.^2)',zeta,wn,wn);
                        end
                end
                
                if k==1
                    S = Sk;
                else
                    S = [S,'.*',Sk];
                end
            end
            if ~isempty(Pars)
                args = sprintf('%s, ',Pars{:});
                argF = sprintf('@(%s, s) ',args(1:end-2));
            else
                argF = '@(s) ';
            end
            h = str2func([argF S]);
            
        end
        function p = pars(obj)
            %PARS
            p = qpar();
            k = 1;
            if isa(obj.par1,'qpar'), p(1,1) = obj.par1; k=k+1; end
            if isa(obj.par1,'qpar'), p(k,1) = obj.par2; end
        end
    end
    
    %methods(Hidden=true)
    methods
        T = rffel(obj,w,dist)
        T = rffpz(obj,w,pzf,dist)
        T = rffcpz(obj,w,pzf,dist)
    end
    
    methods(Static)
        T = rffmul(t1,t2,dist)
        T = rffutil1(w,phi,zmin,zmax,wmin,wmax,form,pzf,kase)
        Tnew = rffutil3(Tleft,T,Tright,edge,dist)
    end
end

