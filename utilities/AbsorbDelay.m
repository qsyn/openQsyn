function [num_d,den_d]=AbsorbDelay(num,den,h,Ts)
%ABSORBDELAYS Absorbs delays into discretized system
%       AbsorbDelay(num,den,h,Ts) recives a certain
%       num/den vectors, delay time h and sample time 
% 		Ts, and returns the equivalent discrete num/den
%       data with the delay absorbed. Currently only 
%       for certain delays.
	nu=h/Ts;
	[Phi,Gamma,c,d]=Local_c2d(num,den,Ts);
    [m,n] = size(Phi); 
    [m,nb] = size(Gamma);
    A_bar=[Phi,Gamma,zeros(m,nu-nb);zeros(nu,m),diag(ones(1,nu-1),1)];
    B_bar=[zeros(nu-1+n,1);1];
    C_bar=[c,zeros(1,nu)];
	[num_d,den_d]=Local_ss2tf(A_bar,B_bar,C_bar,d,'d');
    % Adjust gain
    K=d+c*inv(Phi)*Gamma;
    num_d=K*num_d;
end