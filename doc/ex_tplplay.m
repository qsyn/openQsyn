%% Exmaple: Playing around with templates
%%
% The qtpl class offers conviniante methods for ajustions and
% modifications. This example covers most of them.
% Let us start by generating two random templates: 

template1 = -90 + 10*rand(200,1) + 1i*rand(200,1);
template2 = -100 - 5i + 1*rand(100,1) + 10i*rand(100,1);

tpl1 = qtpl(1,template1); % assign to frequency w=1 rad/s
tpl2 = qtpl(1,template2); % assign to frequency w=1 rad/s
tpl1.show; 
tpl2.show(gcf,'color',[1 0 0]) 
ngrid;

%%
% Observe that since no uncetrain parameters where specified, these filed 
% assumes the default values - indecies. For example: 
tpl1

%%
% Tempaltes can be manipulated using two principle mechanisms: 
%
% * Complex plane (Nyquist) operations are facilitated by |cpop| 
%
% * Nichols plane operations are facilitated by |shift|

%% Complex Plane Operations
% 
% The qtpl method CPOP( A,B,OPR ) performs the operation described by OPR between
% objects A and B, where at least one of them is a qtpl object. useful
% operations include product (|'*'|) and division (|'/'|).

%%  
% *Example 1*: Multiply the template tpl1 by 2.
tpl3 = cpop(tpl1,2,'*');
tpl3.show(gcf,'color',[0 1 0])

%%
% *Example 2*: Multiply tpl1 and tpl2
tpl4 = tpl1.*tpl2 % equivelent to cpop(tpl1,tpl2,'*');
tpl4.show(gcf,'color',[0.5 0 0.5])

%%
% We see that tpl4 has 20000 points. This is bcause |cpop| treat 
% uncertainties as non-dependent, even if they have the smae name (in this
% case, the defualt name).


