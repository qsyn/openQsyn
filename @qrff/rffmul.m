function [T]=rffmul(t1,t2,dist)

%RFFMUL     multiples two templates in real factored form
%	  
%           [T]=rffmul(t1,t2,dist)
%
%       T: 	column vector with an even number of elements.
%       	The first half the vector T is the low gain template edge, 
%           the second half of the vector T is the high gain template edge. 
%       	Each element is of the form degree + j*dB. Each template edge
%       	is sorted in ascending order with respect to angle.
%       	The sequence of angles is not necessarily contiguous.
%       	Each angle is a multiple of dist [deg].
%       	T is found as the upper and lower edges of the vector
%       	addition (concatenation) of t1 and t2 in the Nichols chart
%       	corresponding to the multiplication of the templates in
%       	the Nyquist chart
%
%       t1,t2: 	column vectors representing templates in real factored form,
%       	of the same structure as T, with the relaxation that for each 
%       	edge the angles  do not have to be sorted, nor do they have
%       	to be contiguous (see rffpz). 
%
%       dist: 	The template edges are computed for angles that are
%       	multiples of dist [deg] which must be  a positive real 
%       	number  such that 360 may be divided by dist without
%       	remainder The default of dist is 1 [deg]. The angular 
%       	distance between neigbouring edge points is a multiple of
%       	dist.
%
%
%       Reference: Gutman, P-O, Baril C, Neumann L: "An algorithm for 
%       computing value sets of uncertain transfer functions in 
%       factored real form." IEEE Transactions on Automatic Control,
%       vol 29, no 6, 1268-1273, June 1994.
%
%       See also RFFPZ, RFFEL, RFFCPZ



% Author: P-O Gutman
% Copyright: El-Op Electro-Optics Industries Ltd, and Novomatic AB 1996
% 126-128 204 134-136 146-148
% Version upgrade: A. & Y. Greenhut


if nargin==0 ,
   disp('  T=rffmul(t1,t2,dist)')
   return;
end;


% Defaults
if ~(exist('t2')==1), t2=[]; end;
if ~(exist('dist')==1), dist=1; end;
if (length(t1)==0)&(length(t2)==0),
   T=[]; return;       %error('at least one vector must non-empty');
end;

if length(t2) == 0, T=t1; return; end;		%rem(x,y)= x - round(x./y).*y
if length(t1) == 0, T=t2; return; end;		% matlab built in rem is wrong
if ( (length(t1) == 1) | (length(t2) == 1) ),	%we have to use x-round(x./y).*y
  x=real(t2) - real(t1); 
  
  if  abs(x - round(x./dist).*dist) < 1e-10, %
    T=t1+t2; return; 
  else
    error('non-compatible angles in t1 and t2, or dist non-compatible')
  end; 
end  

t1=t1(:);		%makes one column vector, reading t1 column wise
t2=t2(:);
n1=length(t1);
n2=length(t2);

if (rem(n1,2)~=0), error('Input template t1 has an odd number of elements!'); end
if (rem(n2,2)~=0), error('Input template t2 has an odd number of elements!'); end

t1min=t1(1:n1/2);		% low gain edge
t1max=t1(n1/2+1:n1);		% high gain edge
t2min=t2(1:n2/2);		% low gain edge
t2max=t2(n2/2+1:n2);		% high gain edge

% sort wrt the angles  
[dummy,ix] = sort(real(t1min));
t1min = t1min(ix);
[dummy,ix] = sort(real(t1max));
t1max = t1max(ix);
[dummy,ix] = sort(real(t2min));
t2min = t2min(ix);
[dummy,ix] = sort(real(t2max));
t2max = t2max(ix);

% check if the angles are equal in t1min and t1max, resp t2min and t2max
if (max(abs(real(t1min)-real(t1max)))) > 1e-10,
  error('t1 low gain and high gain edges do not have the same angles')
end  
if (max(abs(real(t2min)-real(t2max)))) > 1e-10,
  error('t2 low gain and high gain edges do not have the same angles')
end 

% check if angles in t1 and t2 are multples of dist (assuming that the angles
% are not multiples of some multiple of dist)
x = real(t1min);
y = real(t2min);
 
if  (max(abs(x - round(x./dist).*dist)) > 1e-10) | (max(abs(y - round(y./dist).*dist)) > 1e-10),
  error(' angles in t1 or t2 are not multiples of dist'); % 
end  

% if one of the templates has only one element in each edge
if (( n1 == 2) | (n2 == 2) ),
  T = [t1min + t2min; t1max + t2max]; return; 
end  

% in each of t1 and t2, check angular distances,
% check contiguity, and complement missing angles with +Inf*j at  low gain edge,
% and -Inf*j at high gain edge

if (~isempty(find(diff(real(t1min))<(dist/2)))) | (~isempty(find(diff(real(t1min)) < (dist/2)))),
  error('multiple occurances of the same angle in t1 or t2, respectively');
end

E=1e15;

if isempty(ix) %
   T=[];			%
else				% end is at bottom
   
while (find(diff(real(t1min)) > ~isempty((3*dist/2)))),
   ix = find(diff(real(t1min)) > (3*dist/2));
   m1 = length(t1min);
   
   if isempty(ix)			% 
      break					%
   end						%
   
   ang = [(real(t1min(ix(1)))+dist):dist:(real(t1min(ix(1)+1))-dist)]';
   t1min = [t1min(1:ix(1)) ; (ang + j*E*ones(size(ang))) ; t1min( (ix(1)+1):m1 ) ];
   t1max = [t1max(1:ix(1)) ; (ang - j*E*ones(size(ang))) ; t1max((ix(1)+1):m1) ];
end
while (find(diff(real(t2min)) > ~isempty((3*dist/2)))),
  ix = find(diff(real(t2min)) > (3*dist/2));
  m2 = length(t2min);
  
  if isempty(ix)			%
     break					%
  end							%
  
  ang = [(real(t2min(ix(1)))+dist):dist:(real(t2min(ix(1)+1))-dist)]';
  t2min = [t2min(1:ix(1)) ; (ang + j*E*ones(size(ang))) ; t2min((ix(1)+1):m2) ];
  t2max = [t2max(1:ix(1)) ; (ang - j*E*ones(size(ang))) ; t2max((ix(1)+1):m2) ];
end
 
% concatenate low gain edges. Anti-Diagonals have the same angle
m1 = length(t1min);	
m2 = length(t2min);

Tmin = t1min*ones(1,m2) + ones(m1,1)*(t2min.'); 	%Kronecker sum 
Tmin = flipud(Tmin);		%Now diagonals have the same angle
					

Ttemp=j*E*ones(min(size(Tmin)), (m1+m2-1));
for k = (-m1+1):(m2-1),
  v = diag(Tmin,k);
  Ttemp(1:length(v),k+m1)=v;
end
Tlo = ( real(Ttemp(1,:)) + j*min(imag(Ttemp)) ).';
  
	%Tlo = zeros(m1+m2-1,1);	%15 % slower alternative to compute Tlo
	%for k = (-m1+1):(m2-1),
	%  v = diag(Tmin,k);
	%  Tlo(k+m1) = real(v(1)) + j*min(imag(v));
	%end

ix = find(imag(Tlo)>1e10);	%eliminate elements with infinite gain
if ~isempty(ix), Tlo(ix)=[]; end

   
% concatenate high gain edges. Anti-Diagonals have the same angle

%m1 = length(t1max);	
%m2 = length(t2max);
Tmax = t1max*ones(1,m2) + ones(m1,1)*(t2max.'); 	%Kronecker sum   
Tmax = flipud(Tmax);

Ttemp= -j*E*ones(min(size(Tmax)), (m1+m2-1));	%faster alternative to compute Thi
for k = (-m1+1):(m2-1),
  v = diag(Tmax,k);
  Ttemp(1:length(v),k+m1)=v;
end
Thi = ( real(Ttemp(1,:)) + j*max(imag(Ttemp)) ).';
 
	%Thi = zeros(m1+m2-1,1);		%slower alternative
	%for k = (-m1+1):(m2-1),
	%  v = diag(Tmax,k);
	%  Thi(k+m1) = real(v(1)) + j*max(imag(v));
	%end

ix = find(imag(Thi)<-1e10);	%eliminate elements with -infinite gain
if ~isempty(ix), Thi(ix)=[]; end

T=[Tlo(:);Thi(:)];
end
%break

%test
%t1 = [(3+j) (6-j) (9-10*j) (12+2*j) (15+0*j) (3+7*j) (6+8*j) (9+10*j) (12+12*j) (15+10*j)].';
%t2 = [(12-2*j) (15-7*j) (18-5*j) (21-4*j) (12+8*j) (15+7*j) (18+11*j) (21+7*j)].';
%dist=3;
