function [numOfconj, conjPairs] = findConjPlace(p,z,flag)
k = 0;
m = 1;
A = [];
switch flag
    case 'pole' 
        for i = 1 : length(p)
            I = imag(p(i));
            for j = 1 : length(p)
                J = imag(p(j));
                if I == - J && I ~=0 
                    A(m,1) = i;
                    A(m,2) = j;
                    k = k + 1;
                    m =  m+1;
                end
            end
        end
    case 'zero'
        for i = 1 : length(z)
            I = imag(z(i));
            for j = 1 : length(z)
                J = imag(z(j));
                if I == - J && I~= 0
                    A(m,1) = i;
                    A(m,2) = j;
                    k = k + 1;
                    m =  m+1;
                end
            end
        end
end
numOfconj = k/2;
conjPairs = A.';
end