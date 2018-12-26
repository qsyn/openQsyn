function [varargout] = bode(obj, w_op, varargin)
%bodetpl  draws a Bode plot from template file
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
%       use PARAMETER='color' with VALUE a color array in RGB format. 
%       use PARAMETER='case' with VALUE equals 'nom', 'unc' (def), or 'all' 
%       use PARAMETER='show' with VALUE equals 'mag', 'phase', or 'both' (def)
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
% 
%   option          string 
%                     'mag'   plot magnitude only
%                     'phase' plot phase only 
%                     empty   plot mag+phase
%

% -------------------------- manage inputs --------------------------------
if nargin==0
   disp('  [phandle]=bodetpl(tplf,w_op,phandle,col) ')
   disp('  [mag,phase,w]=bodetpl(tplf,w_op)')
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
    error('case must be either ''nom'', ''unc'', or ''all''.');
end

if ~(strcmp(opt,'both') ||  strcmp(opt,'mag') || strcmp(opt,'phase'))
    error('show must be either ''mag'', ''unc'', or ''all''.');
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

t_nom=zeros(1,length(w_op));
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
    t_nom(k) = obj(k).template(1);
end

if nargout==3 % get values and do not plot!    
    if strcmp(CASE,'all'), TPL=[t_nom ; TPL]; end
    if strcmp(CASE,'nom'), TPL=t_nom; end
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


% magnitude
if strcmp(opt,'both'), ha(1)=subplot(2,1,1); end
if ~strcmp(opt,'phase')
    if strcmp(CASE,'all') || strcmp(CASE,'unc')
        for k=1:length(tpl)
            semilogx(w_op,imag(TPL(k,:)),'Color',col(k,:),'Tag',num2str(k)); hold on
        end  
    end
    if strcmp(CASE,'all') || strcmp(CASE,'nom')
        [w_nom,t_nom]=gettpl(tplf,'nom');
        semilogx(w_nom,imag(t_nom),'--k','linewidth',2);
    end
    xlim([w_op(1) w_op(end)])
    ylabel('Mag [db]')
    
end

% phase
if strcmp(opt,'both'), ha(2)=subplot(2,1,2); end
if ~strcmp(opt,'mag')
    if strcmp(CASE,'all') || strcmp(CASE,'unc')
        for k=1:length(tpl)
            phase=real(TPL(k,:));
            phase=unwrap(phase*pi/180)*180/pi;
            if phase(1)>5, phase=phase-360; end
            semilogx(w_op,phase,'Color',col(k,:),'Tag',num2str(k)); hold on
        end
    end
    if strcmp(CASE,'all') || strcmp(CASE,'nom')
        phase=real(t_nom);
        phase=unwrap(phase*pi/180)*180/pi;
        if phase(1)>5, phase=phase-360; end
        semilogx(w_nom,real(t_nom),'--k','linewidth',2);
    end
    xlim([w_op(1) w_op(end)])
    ylabel('Phase [deg]')
end

xlabel('Frequency [rad/s]')

% if ~exist('ha'), ha=gca; end
% title(ha(1),getfrom([tplf,'.tpl'],'IO'))

if nargout==1
    varargout{1}=phandle;
end


dcm = datacursormode(gcf);
set(dcm,'updatefcn',@datatipfunc)
end

% custom data cursor tip 
function output_txt = datatipfunc(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
dts = get(event_obj.Target,'Tag');
output_txt = {['Freq: ',num2str(pos(1),4), '[rad/s]'],...
              ['Mag: ',num2str(pos(2),4),' [db]'],...
              ['par set: ',dts]
              };
end

