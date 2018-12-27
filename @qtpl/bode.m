function [varargout] = bode(obj, w_op, varargin)
%bodetpl  draws a Bode plot from template array
%
%   bode(QTPL)   draws the Bode plot of the template in tplf for
%   frequencies given in the file
%
%   bode(QTPL,W)   draws the Bode plot for frequencies given in W
%
%   bode(QTPL,W,FHANDLE)   draws the Bode plot in figure with handle
%   given as FHANDLE
%
%   bode(QTPL,W,...,PARAMETER,VALUE)   parameter/value pairs to 
%   specify additional properties:
%       'color':    VALUE is a color array in RGB format. 
%       'case'      VALUE either 'nom' | 'all' (def)
%       'show'      VALUE either 'mag' | 'phase' | 'magphase' (def)
%   
%   h=bode(...)   returns the handle to the figure 
%
%   [mag,phase,w]=bode(tplf,w_op)   returns matrices for magnitude(db), 
%   phase(deg)and a vector w. Each column in mag,phase corresponds a different 
%   parameter set. No plot is generated.
%
% Output
%
%   phandle         handle to the figure (input+output).
% 
%   mag             Magnitudes in db, n*k matrix
%   phase           Phase in degrees, n*k matrix
%   w               frequencies in rad/s, column vector of length n  
%
% Input
%
%   tplf            Name of template file, from which templates are to be
%                   displayed.
%
%   w_op            Vector of frequencies [rad/s], for which templates are 
%                   to be displayed. A frequency without template is
%                   ignored. If w_op==[], all templates in tplf are
%                   presented.
%   
%   col             Color array in RGB format. if empty, colors are chosen
%                   automatically.

% -------------------------- manage inputs --------------------------------
if nargin==0
   disp('  [phandle]=bode(qtpl,w_op,...) ')
   disp('  [mag,phase,w]=bode(qtpl,w_op)')
   return
end

% default options 
opt='both';
CASE='unc';
col = distinguishable_colors(8);

if nargin>2
    k=1;
    if ishandle(varargin{1}), phandle=varargin{1}; k=2; end
    if isempty(varargin{1}), k=2; end
    
    while k<=nargin-2
        if ischar(varargin{k})
            if nargin<(k+1), error('incorrect parameter/value pairing!'); end
            PARAMETER=varargin{k};
            switch PARAMETER
                case 'color', col=varargin{k+1};
                case 'show', opt=varargin{k+1};
                case 'case', CASE=varargin{k+1};
                otherwise, error('unknown parameter %s',PARAMETER);
            end
        else
            error('incorrect parameter/value pairing!');
        end
        k=k+2;
    end
end

if ~(exist('w_op')==1),  w_op=[];  end
if ~(exist('phandle')==1),  phandle=[];  end

if ~(strcmp(CASE,'nom') ||  strcmp(CASE,'unc') || strcmp(CASE,'all'))
    error('case must be either ''nom'' or ''all''.');
end

if ~(strcmp(opt,'both') ||  strcmp(opt,'mag') || strcmp(opt,'phase'))
    error('show must be either ''mag'', ''phase'', or ''both''.');
end

if isnumeric(col)
    if size(col,2)~=3
        error('''color'' musy be an RGB array');
    end
else
    error('''color'' musy be an RGB array');
end

ncol=size(col,1);
if isempty(w_op), w_op=[obj.frequency]; end
% -------------------------------------------------------------------------

% insert all tpl into a single array
for k=1:length(w_op)
    tpl = obj(k).template(2:end);
    if isempty(tpl), continue; end
    if k==1
        ntpl=length(tpl);
        TPL=nan(ntpl,length(w_op));
    else
        if length(tpl)>ntpl
            tpl=tpl(1:ntpl);
        elseif length(tpl)<ntpl
            ntpl=length(tpl);
            TPL=TPL(1:ntpl,:);
        end
    end
    TPL(:,k)=tpl; % [w1 w2 ... wn]
end

if nargout==3 % get values and do not plot!    
    if strcmp(CASE,'nom'), TPL=TPL(1,:); end
    varargout{1}=imag(TPL).';
    varargout{2}=real(TPL).';
    varargout{3}=w_op;
    return
end 
       
if isempty(phandle)
    phandle=figure('Name','Bode from template','NumberTitle','off');
else
    figure(phandle);
    hold on
end

if length(tpl)>ncol
    col=repmat(col,ceil(length(tpl)/ncol),1);
end


obj.bodeplotter(TPL(2:end,:),w_op,opt,col); % plot all points in template
obj.bodeplotter(TPL(1,:),w_op,opt,[0 0 0]); % plot the nominal case 

if nargout==1
    varargout{1}=phandle;
end

end