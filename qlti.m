classdef qlti < handle
    %QFTTF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        num
        den
        pars
    end
    
    methods
        function G = qlti(varargin)
            
            if nargin==1
                G.num = varargin{1}.num;
                G.den = varargin{2}.den;
                G.pars = varargin{3}.pars;
            end
            if nargin==2
                G.num={};
                G.den={};
                G.pars={};
                G.qftlti2(varargin{1},varargin{2}); 
            end

        end
        function qlti2(Obj,A,B)
            if isa(B,'lti') % Const * H(s)
                H = tf(B);
                [b,a]= tfdata(H);
                b = b{1}; a=a{1};
                Obj.pars = A.pars;
                for k=1:length(b)
                    if b(k)==0
                        Obj.num{k} = 0;
                    elseif b(k)==1
                        Obj.num{k} = sprintf('%s',A.expression);
                    elseif b(k)==-1
                        Obj.num{k} = sprintf('-(%s)',A.expression);
                    else
                        Obj.num{k} = sprintf('%g*%s',b(k),A.expression);
                    end
                end
                Obj.den = num2cell(a);
            end
        end
        function G = plus(A,B)
            G = qftlti;
            if isa(B,'qftlti') % A(s) + B(s)
                G.pars = unique([A.pars; B.pars]);
                na = length(A.num);
                nb = length(B.num);
                G.den=cell(max(na,nb),1);
                A.num = fliplr(A.num);
                A.den = fliplr(A.den);
                B.num = fliplr(B.num);
                B.den = fliplr(B.den);
                for k=1:max(na,nb)
                    if k>na
                        G.num{k} = B.num{k};
                    elseif k>nb
                        G.num{k} = A.num{k};
                    else
                        if A.num{k}==0
                            G.num{k} = B.num{k};
                        elseif B.num{k}==0
                            G.num{k} = B.num{k};
                        else
                            G.num{k} = sprintf('(%s)+(%s)',A.num{k},B.num{k});
                        end
                    end
                end
                G.num = fliplr(G.num);
                G.den = fliplr(G.den);
            end
        end
        %function r = frd(obj)
            
    end
    
end

