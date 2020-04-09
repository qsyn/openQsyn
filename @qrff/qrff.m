classdef qrff
    %QRFF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type
        pars
    end
    
    methods
        function obj = qrff(type,par1,par2)
            %QRFF Construct an instance of the qrff class
            %   Detailed explanation goes here
            
            validStrings = {'gain','delay','dc','hf','uns','int'};
            
            if nargin == 3
                pars = [par1 ; par2];
            elseif nargin ==2
                pars = par1;
            else
                error('Not enough input arguments.')
            end
            if ischar(type)
                validatestring(type,validStrings);
                obj.type = type;
            else
                error('type must be a string')
            end
            if isa(pars,'qpar')
                obj.pars = pars;
            else
                error('pars must be qpar elements')
            end
            
            
            
            
        end
        
        function h = qrff2func(obj)
            %QRFF2FUNC return an handle to a function object f@(p1,p2,...pn,s)
            %with p1,...,pn corresponding to uncertain parameters.
            
            Pars = [];
            for k = 1:length(obj)
                switch obj(k).type
                    case 'gain'
                        Sk = obj(k).pars.name;
                    case 'delay'
                        Sk = sprintf('exp(-s*%s)',obj(k).pars.name);
                    case 'dc'
                        if length(obj(k).pars)==1
                            Sk = sprintf('(1+s/%s)',obj(k).pars.name);
                        else
                            SdW = sprintf('s./%s',obj(k).pars(1).name);
                            Sk = sprintf('(1 + 2.*%s.*(%s) + (%s).^2)',obj(k).pars(2).name,SdW,SdW);
                        end
                    case 'hf'
                        if length(obj(k).pars)==1
                            Sk = sprintf('(s+%s)',obj(k).pars.name);
                        else
                            zeta = obj(k).pars(2).name;
                            wn = obj(k).pars(1).name;
                            Sk = sprintf('(s.^2 + 2.*%s.*%s + %s.^2)',zeta,wn,wn);
                        end
                end
                Pars = vertcat(Pars,obj(k).pars);
                if k==1
                    S = Sk;
                else
                    S = [S,'.*',Sk];
                end
            end
            args = sprintf('%s, ',Pars(:).name);
            argF = sprintf('@(%s, s) ',args(1:end-2));
            h = str2func([argF S]);
            
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

