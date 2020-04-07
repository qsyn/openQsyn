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
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    %methods(Hidden=true)
    methods
        T = rffel(element,a,w,dist)
        T = rffpz(obj,w,pzf,dist) 
        T = rffmul(t1,t2,dist)
        T = rffutil1(w,phi,zmin,zmax,wmin,wmax,form,pzf,kase)
        Tnew = rffutil3(Tleft,T,Tright,edge,dist)
    end
end

