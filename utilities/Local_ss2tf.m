function [num,den]=Local_ss2tf(a,b,c,d,type)
%LOCAL_SS2TF Converts state-space model to num/den.
%       Local_ss2tf(a,b,c,d,type) converts a 
%       state-space model (A,B,C,D) to transfer function
% 		numerator and denominator, 'type' is for discrete
% 		or continuous systems.
%       Problem with gain correction! -- 09/12/29
if isempty(a)
    num=d;
    den=1;
    return
end
	P=eig(a);
	Z=sszero(a,b,c,d);
	switch type
		case 'c'
			if find(P==0,1,'first') %calculate 
				K=c*b;
			else 
				K=d-c*inv(a)*b;
			end
		case 'd'
			if find(P==0,1,'first') %calculate 
				K=c*b;
			else 
				K=d+c*inv(eye(size(a,1))-a)*b;
			end

    end
    % The gain doesn't work, manual fix?
    K=1;
	% Transform back to num/den
	num=K*real(poly(Z)); %maybe need to fix orders
	den=real(poly(P));
end