function [T]=rffutil1(w,phi,zmin,zmax,wmin,wmax,form,pzf,kase)

%RFFUTIL1   utility function for RFFCPZ to compute template edge points
%           according to given edge case for a complex pole or zero pair.
%
%       [T]=rffutil1(w,phi,zmin,zmax,wmin,wmax,form,pzf,case)	
%	
%		
%       T: 	vector of edge template points of the same
%           dimension as phi. Each element is of the form 
%           degree + j*dB, with degree in the interval (-180,180] deg
%	
%       w:  frequency [rad/s], for which the template is computed.
%
%       phi:    column vector with phase values [deg] for which the 
%           edge points are computed.  
%	
%	
%       zmin:	minimum relative damping
%       zmax:	minimum relative damping
%       wmin:	minimum natural frequency [rad/s]
%       wmax: 	maximum natural frequency [rad/s] 
%
%	
%       form:   = 'dc' or 'hf', dc/hf flag
%       pzf: 	= 'p' or 'z',  pole/zero flag
%	
%       case: 	= 1, 2, 3, or 4, border segment case.  
%
%
%       See also RFFCPZ.
%


% Author: P-O Gutman
% Copyright: El-Op Electro-Optics Industries Ltd, and Novomatic AB 1996
% Version upgrade: A. & Y. Greenhut


  if kase == 1
     zx = zmin;
  elseif kase == 3
     zx = zmax;
  elseif kase == 2
     wx = wmax;
  elseif kase == 4
     wx = wmin;
  end 
  
  s = j*w;

if ((kase == 1) | (kase == 3))        % zeta is constant   
  aa=tan(phi*pi/180);   %note that in matlab, tan(pi/2) = 1.6332e+016
  bb=-2*zx*w;  
  cc=-tan(phi*pi/180)*w^2;     

  r1=ones(size(phi))*wmin;	
  r2=ones(size(phi))*wmax;

  if (abs(zx)<eps),     %then bb is small, and phi=0 or 180 deg  
    r1= max(r1); r2=max(r2);    % (see rffcpz) and hence aa=0
  else  
    index3=find(abs(aa)>=1/eps);  %degenerate second order solutions
    index2=find((abs(aa) > eps) & (abs(aa)<1/eps));  % second order solutions
    index1=find((abs(aa) <= eps));  %first order solutions
    r1(index3)=w*ones(size(r1(index3))); %r1(index3)=-w is unnecessary	
    r1(index2)=(-bb+sqrt(bb.^2-4*aa(index2).*cc(index2)))./(2*aa(index2)); 
    r2(index2)=(-bb-sqrt(bb.^2-4*aa(index2).*cc(index2)))./(2*aa(index2));
    r1(index1)=-cc(index1)/bb;
    r2=r2(index2);  		%only index2 contributes to a second sln
  end  
  

  wn=[r1 r2];
  wn=wn(find( (wn>=(wmin-1e-10))&(wn<=(wmax+1e-10)) ));
  z = ones(size(wn))*zx;
  
  if ~isempty(wn)
    % Identify if 'dc' or 'hf'
    %=========================
    if strcmp(form,'dc'),
       t=((s^2)./(wn.^2)+2*z./wn.*s+1);  
    else   
       t=(s^2+2*z.*wn.*s+wn.^2);
    end;

    if abs(zx)<eps,		%wn contains  2 elements only, and phi
                                % contains 0, 1 or 2 elements (0 or 180)
        T=[];
        if ( ((wmin-w)>=0) ),
          if (find(abs(phi)<1e-10)),
            T =  [0     0] + j*imag(c2n(t)); 
          end;
        elseif  ( ((wmax-w)<=0) ),
          if (  find(abs(phi-180)<1e-10) ),  
            T =  [180 180] + j*imag(c2n(t));
          end;
        else                      %wmin<0<wmax
          if (find(abs(phi)<1e-10)),
            T =  [0] + j*imag(c2n(t(2)));
          end;
          if ( find(abs(phi-180)<1e-10)),
            T = [ ([180]+j*imag(c2n(t(1))))    T ];
          end;
        end;  
          
    else  
       T = c2n(t,0);
    end;
  else
    T=[];
  end
     
     
else                                  % wn is constant (case 2 or 4)
  
  z =[];
  
  %make phase in the interval (-180,180] deg
  i180 = find(phi>180);
  if ~isempty(i180), phi(i180) = phi(i180) - 360;	end
  i180 = find(phi<=-180);
  if ~isempty(i180), phi(i180) = phi(180) + 360;  end

  index = find( abs(abs(phi) - 90) > 1e-10);  %~= +-90  (eps does not work)
  if ~isempty(index)
     phi1 = phi(index);
     z =(wx^2-w^2)*tan(phi1*pi/180)/(2*w*wx);
     z=z(find( (z>=(zmin-1e-10))&(z<=(zmax+1e-10)) ));     
  end
  indexvr = find(abs(abs(phi) - 90) < 1e-10);   %= +-90 deg
  if ~isempty(indexvr)
     z=[zmin zmax z];    %detta ar inte ratt i fallet internt phi
  end
  
  wn = ones(size(z))*wx;
  
  % Identify if 'dc' or 'hf'
%=========================
  if strcmp(form,'dc'),
     t=((s^2)./(wn.^2)+2*z./wn.*s+1);  
  else   
     t=(s^2+2*z.*wn.*s+wn.^2);
  end;

  if ~isempty(t)
     T = c2n(t,0);
  else
     T=[];
  end
  
  if ~isempty(indexvr)		% +90 or -90 deg: [zmin zmax] was included 
    indexp90 = find(abs(phi-90)<1e-10);		% +90 deg
    indexm90 = find(abs(phi+90)<1e-10); 	% -90 deg
    if ((~isempty(indexp90)) & (isempty(indexm90))),           %+90 deg only
      if (zmin>=0),
        T(1) = 90 + j*imag(T(1));
        T(2) = 90 + j*imag(T(2));
      elseif (zmax<=0),
        T(1:2) = [];
      else			%zmin<0<zmax
        T(1) = []; 
        T(2) = 90 + j*imag(T(2));
      end;
    elseif ((isempty(indexp90)) & (~isempty(indexm90))),    %-90 deg only
      if (zmin>=0),
        T(1:2) = [];
      elseif (zmax<=0),
        T(1) = -90 + j*imag(T(1));
        T(2) = -90 + j*imag(T(2));  
      else			%zmin<0<zmax
        T(1) = -90 + j*imag(T(1));
        T(2) = [];  
      end  
    elseif ((~isempty(indexp90)) & (~isempty(indexm90)))  %+90 and -90 deg
      if (zmin>=0),
        T(1) = 90 + j*imag(T(1));
        T(2) = 90 + j*imag(T(2));
      elseif (zmax<=0),
        T(1) = -90 + j*imag(T(1));
        T(2) = -90 + j*imag(T(2));  
      else			%zmin<0<zmax
        T(1) = -90 + j*imag(T(1));
        T(2) = +90 + j*imag(T(2));
      end      
    end
  end  
end  
  
%make phase in the interval (-180,180] deg
i180 = find(real(T)>180);
if ~isempty(i180), T(i180) = T(i180) - 360; end
i180 = find(real(T)<=-180);
if ~isempty(i180), T(i180) = T(i180)  + 360; end



% In case   the template segment includes the origin, which then should
% be the low gain edge (cz case), eliminate the origin, since the low gain
% edge will be treated otherwise in RFFCPZ in this case.
% =======================================================

i0 = find(abs(real(T)) < 1e-10);		% 0 deg
if length(T(i0))>1,			%at most 2 elements only if axis is edge
  if (imag(T(i0(1)))<imag(T(i0(2)))) & (imag(T(i0(1)))<imag(c2n(100*eps))),
    T(i0(1))=[];
  elseif   (imag(T(i0(2)))<imag(T(i0(1)))) & (imag(T(i0(2)))<imag(c2n(100*eps))),
    T(i0(2))=[];
  end
end
i0 = find(abs(real(T)-90) < 1e-10);	% 90 deg
if length(T(i0))>1,			%at most 2 elements only if axis is edge
  if (imag(T(i0(1)))<imag(T(i0(2)))) & (imag(T(i0(1)))<imag(c2n(100*eps))),
    T(i0(1))=[];
  elseif   (imag(T(i0(2)))<imag(T(i0(1)))) & (imag(T(i0(2)))<imag(c2n(100*eps))),
    T(i0(2))=[];
  end
end
i0 = find(abs(real(T)-180) < 1e-10);	% 180
if length(T(i0))>1,			%at most 2 elements only if axis is edge
  if (imag(T(i0(1)))<imag(T(i0(2)))) & (imag(T(i0(1)))<imag(c2n(100*eps))),
    T(i0(1))=[];
  elseif   (imag(T(i0(2)))<imag(T(i0(1)))) & (imag(T(i0(2)))<imag(c2n(100*eps))),
    T(i0(2))=[];
  end
end
i0 = find(abs(real(T)+90) < 1e-10);	% -90 deg
if length(T(i0))>1,			%at most 2 elements only if axis is edge
  if (imag(T(i0(1)))<imag(T(i0(2)))) & (imag(T(i0(1)))<imag(c2n(100*eps))),
    T(i0(1))=[];
  elseif   (imag(T(i0(2)))<imag(T(i0(1)))) & (imag(T(i0(2)))<imag(c2n(100*eps))),
    T(i0(2))=[];
  end
end
  


   



