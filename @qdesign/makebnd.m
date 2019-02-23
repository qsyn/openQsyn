function [bound] = makebnd(tpl,specfunc,spec,gphase,gmag)
%MAKEBND calculates a bound 
%   
%   bound = makebnd( tpl,specfunc,spec,gphase,gmag )
%   
%   obj         QDESIGN object (scalar)
%   specfunc    handle to specification criterion function
%   spec        specification bound [upper lower] in dB
%   gphase      grid points for phase
%   gmag        grid points for magnitude
%
% Input check is NOT performaed because this function is not intended to be
% called by the (unexperienced) user

[GPdegM,GPdBM]=meshgrid(gphase,gmag);
GP=GPdegM(:).'+1i*GPdBM(:).';                % A long row vector

GPlength=numel(GP);
tlength=length(tpl);
bnd1=zeros(1,GPlength);

% Is this needed??
N=5000; 
GPsteps=round(linspace(1,GPlength,ceil(min(GPlength*tlength/N,GPlength)))); 
if length(GPsteps)==1, GPsteps = [1 GPlength]; end % enligt Mattias 960221
% --

% who cares??
%if length(GPsteps)>100; 
%	disp(['MAKEBND: Large problem, divided into ',num2str(length(GPsteps)-1),' subproblems']);
%end;

for k=1:(length(GPsteps)-1)
    if k==1  % first interval
        ind=1:GPsteps(2);
    else
        ind=(GPsteps(k)+1):GPsteps(k+1);
    end
    GPtmp=ones(tlength,1)*GP(ind);
    TPLtmp=tpl*ones(1,length(ind));
    %if strcmp(kase,'mimo2')	  %mattias 961102,
    %	PARtmp=par*ones(1,length(ind));
	%else
	%	PARtmp=par;
	%end;
	%if nargin<7,  % no optional parameters to the specification function
	%	bnd1(ind)=feval(specfun,tpl_nom,TPLtmp,GPtmp,spec);
	%else          % optional specification parameters used
	%	bnd1(ind)=feval(specfun,tpl_nom,TPLtmp,GPtmp,spec,par_nom,PARtmp);
	%end;
    bnd1(ind) = specfunc(tpl(1),TPLtmp,GPtmp,spec);
end
bnd=zeros(size(GPdegM));
bnd(:)=bnd1;  % Get back to matrix form suitable for countourc
bound=extract(contourc(gphase,gmag,bnd,[0 5]));

end

function [bnd] = extract(cont)
%EXTRACT    extracts bounds from a contour given by the Matlab command CONTOURC
%
%           [bnd]=extract(cont);
%
%           Subroutine to CBND and BNDUPD
%
%
%           Output:
%
%           bnd     extracted bound vector, with elements in Nichols form
%
%           Input:
%
%           cont    contour vector produced by the Matlab command CONTOURC

% Author: M Nordin
% Version Upgrade: A. & Y. Greenhut
bnd=[];

if isempty(cont)
	bnd = NaN;
	%disp('Warning: No bounds found. Check the search area.');
elseif cont(1,1)~=0
	bnd = NaN;
	%disp('Warning: No bounds found. Check the search area.');
else
	pairs=cont(2,1);
	bnd=cont(1,2:pairs+1)+1i*cont(2,2:pairs+1);  %peo 960411 !
		counter=pairs+2;
	if length(cont)>=counter
		level=cont(1,counter);
		pairs=cont(2,counter);
	else
		level=1;
    end
	while level==0    %obs, foljande rad andrad av Mattias 960411 men fel, sa peo andrade tillbaka
		bnd=[bnd,NaN,cont(1,(counter+1):(counter+pairs))+1i*cont(2,(counter+1):(counter+pairs))];		
		counter=counter+pairs+1;
		%bndsize=size(bnd);
		if length(cont)>=counter
			level=cont(1,counter);
			pairs=cont(2,counter);
		else
			level=1;
        end	
	end % while
end %if

end


