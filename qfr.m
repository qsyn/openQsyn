classdef qfr
    %QFR frequency response function for qft control
    
    properties
        frequency   % frequency in rad/s
        response    % response in Nichols form (deg+j*db)
    end
    
    methods
        function obj = qfr(varargin)
            %QFR Construct an instance of the qfr class
            %
            % Usage:
            %
            %   F = QFR(G)   constructs a qfr object F from the frd object G
            %
            %   F = QFR(G,w)   constructs a qfr object F from the LTI object
            %   G and the freqeucny vector w
            %
            %   F = QFR(response,w)   constructs a qfr object F with given response
            %   and frequency vector
            %
            %   F = QFR(num,den,w)  constructs a qfr object F from a transfer
            %   function given by num and den, and the frequency vector w
            
            % some complex input handling...
            switch nargin
                case 1 % input must be frd
                    if isa(varargin{1},'frd') && isscalar(varargin{1})
                        w = varargin{1}.Frequency;
                        res = c2n(squeeze(varargin{1}.ResponseData));
                    else
                        error('A single input argument to qfr must be a siso frd')
                    end
                case 2 % input is either lti+freq or a freq+response 
                    if isa(varargin{1},'lti') && isscalar(varargin{1}) ...
                            && isnumeric(varargin{2}) && isvector(varargin{2})
                           res = c2n(squeeze(freqresp( varargin{1},varargin{2})));
                           w = varargin{2};
                    elseif isnumeric(varargin{1}) && isvector(varargin{1}) ...
                             && isnumeric(varargin{2}) && isvector(varargin{2})
                         if length(varargin{1})==length(varargin{1})
                             res = varargin{1};
                             w = varargin{2};
                         else
                             error('qfr(response,frequecny) -- both inputs must be of same size')
                         end
                    else
                        error('qfr(response,frequecny) or qfr(G,frequecny) but soething is wrong')
                    end
                case 3
                    if isnumeric(varargin{1}) && isvector(varargin{1}) ...
                            && isnumeric(varargin{2}) && isvector(varargin{2}) ...
                            && isnumeric(varargin{3}) && isvector(varargin{3})
                        
                        num = polyval(varargin{1},1i*varargin{3});
                        den = polyval(varargin{2},1i*varargin{3});
                        res = c2n(num./den);
                        w = varargin{3};
                    else
                        error('qfr(num,den,freqeucny) but something is wrong')
                    end
                otherwise
                    error('inputs are not good')
            end
            % make sure properties are column vector for consistency
            obj.response = reshape(res,[],1);
            obj.frequency = reshape(w,[],1);
        end
        function varargout = show(obj,varargin)
            %SHOW show response on Nichols chart
            h = plot(real(obj.response),imag(obj.response),varargin{:});
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
        function [] = bodeplot(obj,varargin)
            %BODEPLOT plots a bode plot with additional options
            G = obj.frd();
            opt = bodeoptions(varargin{:});
            bodemag(G,obj.frequency,opt);
        end
        function G = frd(obj)
            %FRD converts qfr to a Control Systems Toolbox FRD object
            res = n2c(obj.response);
            G = frd(res,obj.frequency);    
        end   
        function G = series(A,B)
           %SERIES conection of QFR object with another QFR or LTI object
           w = A.frequency;
           switch class(B)
               case 'qfr'
                   if all(w == B.frequency)
                       G = qfr(A.response+B.response,A.frequency);
                   else
                       error(['series connection of QFR object require that both ',...
                           'objects have identical frequency vector']); 
                   end
               case {'tf','zpk','ss','frd'}
                   Bfr = squeeze(freqresp(B,w)).';
                   Bnic = c2n(Bfr.',-1);
                   G =  qfr(A.response+Bnic,w);
               case 'double'
                   %if ~isscalar(B), error('a numeric value must be a scalar'); end
                   G = qfr(A.response+B,w);
               otherwise
                   error('seconf input must be an QFR object, LTI, or a double')
           end
        end
        function t = freqresp(obj,w)
            %FREQRESP returns the frequency resposne at selected frequency
            %If w is not in obj.frequency, FREQRES interpolates
                        
            t = interp1(obj.frequency,obj.response,w);
            
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

