% To test the blackBox option we repeat the plant from the basic siso
% exmaple. It is found in examplePlant.m

f = @examplePlant

% uncertain parameters
k=qpar('k',2,2,5,8);
a=qpar('a',3,1,3,8);
z=qpar('z',0.6,0.3,0.6,8);
wn=qpar('wn',4,4,8,8);

% design frequencies
w = [0.2 0.5 1 2 5 10 20 50];

% make a qpar array as a column vector (horizontal concatenation produces a qpoly object!)
pars=[k ; a ; z ; wn];

Pblack = qplant(f,pars)


Pblack.cnom(logspace(-2,2,200))