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
            [Lia,Locb] = ismember(obj.poles,obj.zeros);
            newobj = obj;
            newobj.poles(Lia) = [];
            newobj.zeros(Locb(Locb>0))= [];
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
        function q = qfr(obj,frequency)
           %QFR converts qctrl object into a qfr object with specified frequency
           res = nicresp(obj,frequency);
           q = qfr(res,frequency);
        end
        function str = print(obj)
            %PRINT prints as string in dc form
            %
            if obj.gain==1
                s1=[];
            else
                s1 = sprintf('%g',obj.gain);
            end
            s2 = [];
            s3 = [];
            if isempty(obj.zeros)
                s2 = '1';
            end
            k=1;
            while k <= length(obj.zeros)
                z = obj.zeros(k);
                if isreal(z)
                    if z==0 
                        s2 = [s2,'s'];
                    elseif z==1
                        s2 = [s2,'(1-s)'];
                    elseif z==-1
                        s2 = [s2,'(s+1)'];
                    elseif z>0
                        s2 = [s2,sprintf('(1-s/%g)',z)];
                    else
                        s2 = [s2,sprintf('(s/%g+1)',-z)];
                    end
                    k=k+1;
                else
                    wn = abs(z);
                    zeta = real(z)/wn;
                    if zeta<0
                        zeta = -zeta;
                        wn = -wn;
                    end
                    s2 = [s2,sprintf('(s^2/%g+%g*s/%g+1)',wn^2,2*zeta,-wn)];
                    k=k+2;
                end
            end
            %if length(s2)>1
            %    s2 = s2(1:end-1);
            %end
            k=1;
            while k <= length(obj.poles)
                p = obj.poles(k);
                if isreal(p)
                    if p==0 
                        s3 = [s3,'s'];
                    elseif p==-1
                        s3 = [s3,'(s+1)'];
                    elseif p==1
                        s3 = [s3,'(1-s)'];
                    elseif p>0
                        s3 = [s3,sprintf('(1-s/%g)',p)];
                    else
                        s3 = [s3,sprintf('(s/%g+1)',-p)];
                    end
                    k=k+1;
                else
                    wn = abs(p);
                    zeta = real(p)/wn;
                    if zeta<0
                        zeta = -zeta;
                        wn = -wn;
                    end
                    s3 = [s3,sprintf('(s^2/%g+%g*s/%g+1)',wn^2,2*zeta,-wn)];
                    k=k+2;
                end
            end
            %if ~isempty(s3)
            %    s3 = s3(1:end-1);
            %end
            [L,I] = max([length(s2) length(s3)]);
            sline = repmat('-',1,L+2);
            pad1 = repmat(' ',1,length(s1));
            pad2 = repmat(' ',1,ceil(abs(length(s2)-length(s3))/2));
            if I==1
                s3 = [pad2,s3,pad2];
            else
                s2 = [pad2,s2,pad2];
            end
            
            str = sprintf('    %s%s\n   %s%s\n   %s %s ',pad1,s2,s1,sline,pad1,s3);
        end
        function B = inv(A)
            %INV inverts a qctrl object
            B = qctrl(A.poles,A.zeros,1/A.gain);
            B.sampleTime = A.sampleTime;
        end
        function C = series(A,B)
            %SERIES connection
            
            if isa(A,'qctrl') && isa(B,'qctrl')
                if A.sampleTime ~= B.sampleTime
                    error('sample time must bbe the same')
                end
                C = qctrl([A.zeros.' B.zeros.'],[A.poles.' B.poles.'],A.gain*B.gain);
                C.sampleTime = A.sampleTime;
            elseif isa(A,'qctrl') && isnumeric(B)
                C = qctrl(A.zeros,A.poles,A.gain*B);
                C.sampleTime = A.sampleTime;
            elseif isnumeric(A) && isa(B,'qctrl')
                C = qctrl(B.zeros,B.poles,B.gain*A);
                C.sampleTime = B.sampleTime;
            else
                error('illigal inputs!')
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
           %k = obj.gain/obj.zeros(end)*obj.poles(end);
           if obj.sampleTime == 0
               T0 = zpk(obj.zeros,obj.poles,1,'displayformat','frequency');
           else
               T0 = zpk(obj.zeros,obj.poles,1,obj.sampleTime,'displayformat','frequency');
           end
           T = obj.gain/dcgain(T0)*T0; % adjust dc gain
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
        obj = lag(Freq,beta)
        obj = notch(wn,wd,zn,zd)
        obj = PID(Phase,Freq,Damping)
        obj = Qpz(a,b,c,d,flag)
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

