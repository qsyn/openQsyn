classdef qfr
    %QFR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frequency
        nic
    end
    
    methods
        function obj = qfr(nresponse,frequecny)
            %QFR Construct an instance of this class
            %   Detailed explanation goes here
            obj.frequency = frequecny;
            obj.nic = nresponse;
        end
        function varargout = show(obj,varargin)
            %NIC plots the on Nichols chart
            h = plot(real(obj.nic),imag(obj.nic),varargin{:});
            if nargout==1, varargout{1}=h; end   
        end
        function [] = nichols(obj)
            %NICHOLS plots a nichols chart, compatible to Control System
            %Toolbox
            G = obj.frd();
            nichols(G,obj.frequency);
        end
        function [] = bode(obj)
            %BODE plots a bode plot
            G = obj.frd();
            bode(G,obj.frequency);
        end
        function G = frd(obj)
            %FRD converts qfr to frd object
            
            response = n2c(obj.nic);
            G = frd(response,obj.frequency);    
        end 
        function G = series(A,B)
           %SERIES conection of QFR object with another QFR or LTI object
           w = A.frequency;
           switch class(B)
               case 'qfr'
                   if all(w == B.frequency)
                       G = qfr(A.nic+B.nic,A.frequency);
                   else
                       error(['series connection of QFR object require that both ',...
                           'objects have identical frequency vector']); 
                   end
               case {'tf','zpk','ss','frd'}
                   Bfr = squeeze(freqresp(B,w)).';
                   Bnic = c2n(Bfr,'unwrap');
                   G =  qfr(A.nic+Bnic,w);
               case 'double'
                   %if ~isscalar(B), error('a numeric value must be a scalar'); end
                   G = qfr(A.nic+B,w);
               otherwise
                   error('seconf input must be an QFR object, LTI, or a double')
           end
        end
        function t = freqresp(obj,w)
            %FREQRESP returns the frequency resposne at selected frequency
            %If w is not in obj.frequency, FREQRES interpolates
                        
            t = interp1(obj.frequency,obj.nic,w);
            
        end
    end
    
    methods(Static)
        function Gqfr =lti2qfr(Glti,Frequency)
            %LTI2QFR converets an LTI object from the Control Systems
            %Toolbox into a QFR object
            if ~isa(Glti,'lti'), error('first input must be an LTI object'); end
            if ~(isnumeric(Frequency) && isvector(Frequency))
                error('second input must be a numeric vector')
            end
            
            Gfr = freqresp(Glti,Frequency);
            Gnic = c2n(Gfr,'unwrap');
            Gqfr = qfr(Gnic,Frequency);
        end
    end
    
end

