function zVec=sszero(A,B,C,D)
%SSZERO Zeros of a state-space model.
%       sszero(A,B,C,D) computes the (transmission) zeros of the
%       state-space model (A,B,C,D).  
maxz=1000;
AA=([A B;C D]);
BB=[eye(length(A)) B*0;0*[C D]];
if length(B(1,:)) ~= length(C(:,1))
 fprintf('\nSSZERO only works for square systems \n')
 return
end
[AA,BB,Q,V,Z]=qz(AA,BB);
BB=diag(BB);
AA=diag(AA);
I=find(maxz*abs(BB) > abs(AA));
zVec=AA(I)./BB(I);
end
