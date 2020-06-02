function str = print(obj)
%PRINT prints QCTRL object as a string in dc form
%
if obj.gain==1
    s1=[];
else
    s1 = sprintf('%g',obj.gain);
end
s2 = [];
s3 = [];
if isempty(obj.zeros)
    s2 = '1';
end
k=1;
%nzeros = sum(obj.zeros==0);     % number of zeros
%k = 1+nzeros;
%if nzeros>1
%    s2 = sprintf('s^%i%',nzeros);
%end
nzeros=0;
while k <= length(obj.zeros)
    z = obj.zeros(k);
    if isreal(z)
        if z==0
            %s2 = [s2,'s'];
            nzeros = nzeros+1;
        elseif z==1
            s2 = [s2,'(1-s)'];
        elseif z==-1
            s2 = [s2,'(s+1)'];
        elseif z>0
            s2 = [s2,sprintf('(1-s/%g)',z)];
        else
            s2 = [s2,sprintf('(s/%g+1)',-z)];
        end
        k=k+1;
    else
        wn = abs(z);
        zeta = real(z)/wn;
        if zeta<0
            zeta = -zeta;
            wn = -wn;
        end
        s2 = [s2,sprintf('(s^2/%g+%g*s/%g+1)',wn^2,2*zeta,-wn)];
        k=k+2;
    end
end
if nzeros>1
    s2 = sprintf('s^%i%s',nzeros,s2);
elseif nzeros>0
    s2 = ['s',s2];
end
%if length(s2)>1
%    s2 = s2(1:end-1);
%end
k=1;
npoles=0;
while k <= length(obj.poles)
    p = obj.poles(k);
    if isreal(p)
        if p==0
            %s3 = [s3,'s'];
            npoles=npoles+1;
        elseif p==-1
            s3 = [s3,'(s+1)'];
        elseif p==1
            s3 = [s3,'(1-s)'];
        elseif p>0
            s3 = [s3,sprintf('(1-s/%g)',p)];
        else
            s3 = [s3,sprintf('(s/%g+1)',-p)];
        end
        k=k+1;
    else
        wn = abs(p);
        zeta = real(p)/wn;
        if zeta<0
            zeta = -zeta;
            wn = -wn;
        end
        s3 = [s3,sprintf('(s^2/%g+%g*s/%g+1)',wn^2,2*zeta,-wn)];
        k=k+2;
    end
end
if npoles>1
    s3 = sprintf('s^%i%s',npoles,s3);
elseif npoles>0
    s3 = ['s ',s3];
end
if isempty(s3)
    str = sprintf('%s %s',s1,s2);
else
    [L,I] = max([length(s2) length(s3)]);
    sline = repmat('-',1,L+2);
    pad1 = repmat(' ',1,length(s1));
    pad2 = repmat(' ',1,ceil(abs(length(s2)-length(s3))/2));
    if I==1
        s3 = [pad2,s3,pad2];
    else
        s2 = [pad2,s2,pad2];
    end
    str = sprintf('    %s%s\n   %s%s\n   %s %s ',pad1,s2,s1,sline,pad1,s3);
end
end

