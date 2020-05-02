%% Computing Horowitz-Sidi bounds
% After all specification are inserted into |qspc| objects, and all templates
% are calculated, we may finally proceed calculating the Horowitz-Sidi bounds.
% This is done using the <matlab:doc('qdesign') |qdesign|> class,
% constructed as
%
%   des = qdesign(P,sepcs);
% 
% where P is a |qplant| object with computed nominal response and
% templates, and |specs| is an array of all |qspc| objects storing the specifications. 
%
% Following, each bound is computed using the command |cbnd|:
%
%    des.cbnd(name)
%
% with |name| a string giving the specification to compute 
%
% To show a computed bound use |des.show(name)|. 
% To plot all dominant bound use |des.show('dom')|.


