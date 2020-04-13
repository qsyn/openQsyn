function [nt]=cltmp(t,dist)

%CLTMP      interpolates points in a sorted and pruned template (subroutine) 
%           [nt]=cltmp(t,dist) 
%
%           closes the distance between two neighbouring points
%           in the template to a desired distance, by adding new template points
%           along the line between them in the Nichols chart. Subroutine to rff.
%	
%
%           Output:
%           
%           nt:     new template vector with elements in Nichols form [deg + j * dB]
%           
%           Inputs:
%
%           t       pruned and sorted template vector with elements in Nichols form 
%                   [deg + j * dB]
%
%           dist     a two elements  row vector 
%                   [ phase distance in degrees , gain distance in dB ]
%                   denoting the maximal  distance (2-norm) between two
%                   neighbouring template points in the new template

%
%           Remark: The template t must be sorted, i.e  the points p(i-1) and p(i+1) 
%   	            are the nearest points to the point p(i) 
%



% Author: B Cohen
% Copyright: El-Op Electro-Optics Industries Ltd, and Novomatic AB 1996
% Version Upgrade: A. & Y. Greenhut
 

if nargin==0 ,
   disp('  [nt]=cltmp(t,dist)')
   return;
end;

t=t(:);
if length(t)==0, nt=[]; return; end;
[n,m]=size(t);
t(n+1,:)=t(1,:);
nt=[];

for i=1:n,   
    dx=real(t(i+1,1)-t(i,1));
    dy=imag(t(i+1,1)-t(i,1));
    if dx ~= 0, ndx=dx/dist(1); else ndx=0; end;
    if dy ~= 0, ndy=dy/dist(2); else ndy=0; end;
    ndiv=max(ceil(abs([ndy+j*ndx])));
    if ndiv > 1 ,     
     xnt=[];
      y=imag([t(i) t(i+1)]);
      x=real([t(i) t(i+1)]);
      if abs(diff(x)) < 1e-10 , 
         yi=linspace(y(1),y(2),ndiv+2);
         yi=yi(2:ndiv+1);
         xi=ones(1,ndiv)*x(1);         
      else
         c=polyfit(x,y,1);
         xi=linspace(x(1),x(2),ndiv+2);
         xi=xi(2:ndiv+1);
         yi=polyval(c,xi);         
      end;         
      t_=xi+yi*j;
      xnt=[xnt t_(:)];      
      nt=[nt ; t(i,:) ; xnt];
    else
         nt=[nt ; t(i,:)];
    end;   
    
end;


