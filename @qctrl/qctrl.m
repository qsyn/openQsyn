classdef qctrl < matlab.mixin.CustomDisplay
    %QCTRL is a class for openQsyn controllers
    %   Use for control design w/o Control System Toolbox
    
    properties
        gain            (1,1)   {mustBeNumeric,mustBeReal}
        poles           (:,1)   {mustBeNumeric}
        zeros           (:,1)   {mustBeNumeric}
        sampleTime      (1,1)   {mustBeNumeric,mustBeReal,mustBeNonnegative}
    end
    
    methods
        function obj = qctrl(varargin)
            %QCTRL Construct an instance of this class
            %   Detailed explanation goes here
            
            if nargin==2 || nargin>4
                error('wrong number of input arguments');
            end
            p=inputParser;
            p.addOptional('input1',[]); % Zeros/Gain/LTI
            p.addOptional('input2',[]); % Poles
            p.addOptional('input3',1);  % Gain
            p.addOptional('input4',0);  % Sample Time
            p.parse(varargin{:});
            
            if nargin ==1 && isa(p.Results.input1,'lti')
                obj.zeros = zero(p.Results.input1);
                obj.poles = pole(p.Results.input1);
                nint = sum(obj.poles==0); % number of integratores
                ndif = sum(obj.zeros==0); % number of deifferentiators
                T0 = tf([1 zeros(1,nint)],[1 zeros(1,ndif)]); % remove integratros/differentiators 
                obj.gain = dcgain(p.Results.input1*T0);
                if isdt(p.Results.input1)
                    obj.sampleTime = p.Results.input1.Ts;
                end
            elseif nargin==1 && isnumeric(varargin{1})
                obj.gain = varargin{1};
            else
                obj.zeros = p.Results.input1;
                obj.poles = p.Results.input2;
                obj.gain = p.Results.input3;
                obj.sampleTime = p.Results.input4;
            end
            
        end
        function newobj = minreal(obj)
            %MINREAL
            %[Lia,Locb] = ismember(obj.poles,obj.zeros);
            newobj = obj;
            k=1;
            while k<=length(newobj.poles)
                [~,I] = ismember(newobj.poles(k),newobj.zeros);
                if I>0
                    newobj.poles(k)=[];
                    newobj.zeros(I)=[];
                else
                    k=k+1;
                end
            end
            %newobj.poles(Lia) = [];
            %newobj.zeros(Locb(Locb>0))= [];
        end
        function [num,den] = tfdata(obj)
            %TFDATA computes numerator and denominator tf data
            Idiff = obj.zeros==0;
            Iint = obj.poles==0;
            num = poly(obj.zeros(~Idiff));
            den = poly(obj.poles(~Iint));
            num = num/num(end)*obj.gain;
            den = den/den(end);
            num = [num zeros(1,sum(Idiff))];
            den = [den zeros(1,sum(Iint))];
        end
        function res = nicresp(obj,frequency)
            %NICRESP returns the frequeny reponse in Nichols format
            
            p = inputParser;
            p.addRequired('frequency',@(x) validateattributes(x,...
                {'numeric'},{'nonnegative','real','finite','increasing','vector'}));
            p.parse(frequency);
            w = p.Results.frequency;
            s = frequency*1i;
            [num,den] = tfdata(obj);
            res = c2n(polyval(num,s)./polyval(den,s),-180);
            
        end
        function nichols(obj,frequency,varargin)
           %NICHOLS plots Nichols chart for qctrl objects
           %
           %Usage: 
           %    NICHOLS(C)   plots Nichols chart for QCTRL object C over
           %    defualt frequency grid
           %
           %    NICHOLS(C,W)   specify the frequency vector
           %
           %    NICHOLS(C,W,LineSpec)   sets the line style, marker symbol, 
           %    and color
           %    
           %    plot(___,Name,Value)   specifies line properties using one 
           %    or more Name,Value pair arguments.
           
           if ~exist('frequency','var'), frequency=[]; end
           if isempty(frequency), frequency=logspace(-2,2,200); end
           
           res = nicresp(obj,frequency);
           fr = qfr(res,frequency);
           show(fr,varargin{:});
                   
        end
        function q = qfr(obj,frequency)
           %QFR converts qctrl object into a qfr object with specified frequency
           res = nicresp(obj,frequency);
           q = qfr(res,frequency);
        end
        function B = inv(A)
            %INV inverts a qctrl object
            B = qctrl(A.poles,A.zeros,1/A.gain);
            B.sampleTime = A.sampleTime;
        end
        function B = uplus(A)
           %UPLUS unitary plus (+)
           B = A;
        end
        function B = uminus(A)
           %UMINUS unitary minus (-)
           B = A;
           B.gain = -B.gain;
        end
        function C = plus(A,B)
            %PLUS addition of qctrl objects, a parallel connection
            C = parallel(A,B);
        end
        function C = minus(A,B)
            %MINUS substruction of qctrl objects, a parallel connection
            C = parallel(A,-B);
        end
        function C = series(A,B)
            %SERIES connection
            
            if isa(A,'qctrl') && isa(B,'qctrl')
                if A.sampleTime ~= B.sampleTime
                    error('sample time must agree')
                end
                C = qctrl([A.zeros.' B.zeros.'],[A.poles.' B.poles.'],A.gain*B.gain);
                C.sampleTime = A.sampleTime;
            elseif isa(A,'qctrl') && isnumeric(B)
                C = qctrl(A.zeros,A.poles,A.gain*B);
                C.sampleTime = A.sampleTime;
            elseif isnumeric(A) && isa(B,'qctrl')
                C = qctrl(B.zeros,B.poles,B.gain*A);
                C.sampleTime = B.sampleTime;
            elseif isa(B,'qsys') || isa(B,'qplant')
                C = series(B,A);
            else
                error('qctrl object may only be connected to objects of class qctrl, qplant, qsys, or numerical scalars')
            end
            
        end
        function C = mtimes(A,B)
            %MTIMES is equivalent to serial connection
            %
            % C = A*B
            C = series(A,B);
            
        end
        function C = mldivide(A,B)
            %MLDIVIDE returns A\B
            C = series(inv(A),B);
        end
        function C = mrdivide(A,B)
            %MRDIVIDE returns A/B
            C = series(A,inv(B));
        end
        function C = mpower(A,b)
           %MPOWER power, ^
           p=inputParser;
           p.addRequired('b',@(x) validateattributes(x,{'numeric'},{'scalar','integer','real'}))
           p.parse(b);
           C = qctrl([],[],1);
           
           for k = 1:abs(p.Results.b)
               C = series(A,C);
           end
           if p.Results.b<0
            C = 1/C;
           end
        end
        function bodeplot(obj,frequency,varargin)
           %BODEPLOT Bode plot of qctrl objects 
           %
           %Usage: BODEPLOT(Q,W)   generates a Bode plot of QCTRL object Q on a
           %given frequency grid W. 
 
           if ~exist('frequency','var')
               frequency = logspace(-2,3,200);
           end
           res = nicresp(obj,frequency);
           QFR = qfr(res,frequency);
           bodeplot(QFR,varargin{:});
        end
        function T = zpk(obj)
           %ZPK converts qctrl to zpk class
           Idiff = obj.zeros==0;
           Iint = obj.poles==0;
           if obj.sampleTime == 0
               T0 = zpk(obj.zeros(~Idiff),obj.poles(~Iint),1,'displayformat','frequency');
           else
               T0 = zpk(obj.zeros(~Idiff),obj.poles(~Iint),1,...
                   obj.sampleTime,'displayformat','frequency');
           end
           % adjust dc gain
           Tadj = obj.gain/dcgain(T0)*T0; 
           T = Tadj*zpk(zeros(1,sum(Idiff)),zeros(1,sum(Iint)),1);
        end
        function T = tf(obj)
            %TF converts qctrl to tf class
            [num,den]=tfdata(obj);
            if obj.sampleTime == 0
                T = tf(num,den);
            else
                T = tf(num,den,obj.sampleTime);
            end
        end
    end
    
    
    methods(Static)
        obj = lead(Phase,Freq,Damping)
    end
    
    methods(Access = protected)
        function displayScalarObject(obj)
            disp(obj.print)
            if obj.sampleTime==0
                fprintf('\nContinuous-time openQsyn QCTRL object\n\n')
            else
                fprintf('\nDiscrete-time openQsyn QCTRL object, sample time: %g\n\n',obj.sampleTime);
            end
        end
    end
    
end

