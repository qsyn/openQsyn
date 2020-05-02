function Clag = Clag(Freq,beta)
wn =  Freq;
s = tf('s');
Clag  = (10*s+wn)/(10*s+wn/beta);
end