function [Phi,Gamma,c,d]=Local_c2d(num,den,Ts)
%LOCAL_C2D Discretizied state-space model.
%       Local_c2d(num,den,Ts) computes the zero-order
%       hold equivalent state-space model (A,B,C,D)
% 		using Van Loan's matrix exponential formula.
	[a,b,c,d] = Local_tf2ss(num,den);
    % Van Loan
    [m,n] = size(a); 
    [m,nb] = size(b); 
    str = expm([[a b]*Ts; zeros(nb,n+nb)]);
    Phi = str(1:n,1:n);
    Gamma = str(1:n,n+1:n+nb);	
end
