function [Tnew]=rffutil3(Tleft,T,Tright,edge,dist)

%RFFUTIL3   utility function for RFFCPZ to sort, clean, and complement edges.
%
%           [Tnew]=rffutil3(Tleft,T,Tright,edge,dist)
%
%		The edge T is first sorted w r t ascending angle.
%		Then the sorted edge, T,  and which has groups of more than one
%		element with the same phase, select, from each group, the highest
%		gain member(s), if edge=='hi', and the lowest gain member if edge=='lo'.
%		Finally, the edge is complemented with one or two new endpoints if these
%		are angularly nearer its true endpoints.
%		
%		Angular elimination may occur for -90, 0, 90, 180 degrees,
%		when an edge sits or "almost" sits on an axis, or at the vertices
%		of elementary edges. Phase rounding in described in RFFCPZ
%
%		Tleft:	Two exact endpoints of left (low angle) segment of T
%		T: 	edge to be corrected, row or column vector in Nichols form
%			(degree + j*dB) 
%               Tright: Two exact endpoints of right (high angle) segment of T
%		edge:	'hi' or 'lo'              
%		dist:	current angular resolution, degrees 
%		Tnew:   sorted, cleaned, and end point complemented edge
%
%		See RFFCPZ


% Author: P-O Gutman
% Copyright: El-Op Electro-Optics Industries Ltd, and Novomatic AB 1996
% Version upgrade: A. & Y. Greenhut

 
  [dummy,ix]=sort(real(T));		%sort
  T=T(ix);
  
  n = length(T);			%clean
  Tnew = T;
  if (n > 0) & strcmp(edge,'hi'),
    while ~isempty(find(abs(diff(real(Tnew)) < (dist/2))))
      ix = find(abs(diff(real(Tnew)) < (dist/2))); 
      if imag(Tnew(ix(1)))<imag(Tnew(ix(1)+1)),
        Tnew(ix(1))=[];
      else
        Tnew(ix(1)+1)=[];
      end
    end
  elseif  (n > 0) & strcmp(edge,'lo'),
    while ~isempty(find(abs(diff(real(Tnew)) < (dist/2))))
      ix = find(abs(diff(real(Tnew)) < (dist/2))); 
      if imag(Tnew(ix(1)))>imag(Tnew(ix(1)+1)),
        Tnew(ix(1))=[];
      else
        Tnew(ix(1)+1)=[];
      end
    end
  end  

  % Complement extreme points with phase rounding
  if (~isempty(Tleft)) | (~isempty(Tright)),
    if ~isempty(Tnew),
      if real(Tleft(1)) < real(Tnew(1)) - dist/2,
        Tnew = [(real(Tnew(1)) - dist + j*imag(Tleft(1))) ; Tnew];
      end     
      n = length(Tnew);
      if length(Tright)==1, Tright(2)=Tright(1); end
      if real(Tright(2)) > real(Tnew(n)) + dist/2,
        Tnew = [Tnew; (real(Tnew(n)) + dist + j*imag(Tright(2)))];
      end;
    else            % true template lies between phase grid lines
      Tedge = [ Tleft Tright];
      Tnewtemp = round(real(Tedge)/dist)*dist + j*imag(Tedge);
      m = length(Tnewtemp);
      if real(Tnewtemp(m)-Tnewtemp(1))<(dist/2),  % same phase rounding
        if edge=='hi', 
          Tnew = real(Tnewtemp(1)) + j*max(imag(Tnewtemp));
        else
          Tnew = real(Tnewtemp(1)) + j*min(imag(Tnewtemp));
        end
      else                                  % adjacent phase rounding
        ileft = find(real(Tnewtemp)==real(Tnewtemp(1))); %left group
        if edge=='hi', 
          Tnew = real(Tnewtemp(1)) + j*max(imag(Tnewtemp(ileft)));
        else
          Tnew = real(Tnewtemp(1)) + j*min(imag(Tnewtemp(ileft)));
        end        
        iright = find(real(Tnewtemp)==real(Tnewtemp(m))); %right group
        if edge=='hi', 
          Tnew = [Tnew; real(Tnewtemp(m)) + j*max(imag(Tnewtemp(iright)))];
        else
          Tnew = [Tnew; real(Tnewtemp(m)) + j*min(imag(Tnewtemp(iright)))];
        end                
      end
    end
  end        
      
  
