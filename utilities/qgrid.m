function[outgrid]=qgrid(n,qmin,qmax);
%
%QGRID      makes a grid out of two vectors 
%
%           [outgrid]=qgrid(n,qmin,qmax);
%
%	First grids rowwise between qmin(i) and qmax(i) with n(i) grid points,
% 	and then produces all vector combinations over i, such that each column 
%   in [outgrid] is a combination.
%
%   Output:
%
%   [outgrid]   matrix where each column is a combination. [outgrid] has
%               as many rows as n has elements. 
%
%   Inputs:
%
% 	n           integer vector. Each element n(i) tells how many grid points
%               there will be in the interval [qmin(i) qmax(i)]
%
%   qmin        vector of the same size as n, denoting the left endpoints
%
%   qmax        vector of the same size as n, denoting the left endpoints
%
%	Example: 
%  	            Q=qgrid([3,2];[1 10],[2 20]), 
%  	                    Q= [  1.0   1.5   2.0    1.0    1.5    2.0
%                            10.0  10.0  10.0   20.0   20.0   20.0 ]
%
%   See also    PGRID, PARGRID


% Author: C Baril, B Cohen
% Copyright: El-Op Electro-Optics Industries Ltd, and Novomatic AB 1996
% Version Upgrade: A. & Y. Greenhut

n=n(:);qmin=qmin(:)';qmax=qmax(:)';  %Make column vectors;
j=0:prod(n)-1;
outgrid=zeros(size(n))*zeros(size(j));
for k=1:length(n);
   outgrid(k,:)=rem(fix(j/prod(n(1:k-1))),n(k))+1;
   if n(k)==1
      outgrid(k,:)=(qmin(k)+qmax(k))*0.5*ones(size(outgrid(k,:)));
   else
      outgrid(k,:)=((n(k)-outgrid(k,:))*qmin(k)+(outgrid(k,:)-1)*qmax(k))*(1/(n(k)-1));
   end;
end;   
   
