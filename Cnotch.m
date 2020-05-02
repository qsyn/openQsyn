function Cnotch = Cnotch(w0n,w0d,zn,zd)
s = tf('s');
Cnotch = ((s^2)/(w0n^2)+2*zn*(s/w0n)+1)/((s^2)/(w0d^2)+2*zd*(s/w0d)+1);
end