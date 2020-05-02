%% Open Qsyn 
%%
% *Open Qsyn* is a modern open source toolbox for _QFT_ control _synthesis_. 
%
% _Quantitative Feedback Theory_ (QFT) is a frequency domain robust control
% design technique, introduced by Isaac Horowitz.
%
% This toolbox is the successor of the original _Qsyn_ toolbox, developed in 
% the 90s. It provides a modern and completely free open source toolbox to 
% aid QFT control synthesis. The development is supported by the author of
% the original _Qsyn_ toolbox, Prof. Per-Olof Gutman, and all reused code is done 
% under his premission.
%
% The project is hosted on GitHub: <https://github.com/qsyn/openQsyn> 

%% 
% Typical QFT design is doen according to the following steps:
% 
% # <matlab:web('doc_modeling.html') Modeling of the uncertain plants>  
% # <matlab:web('doc_tplcomp.html') Computating  templates> -- the set of all plant responses over the Nichols
% chart at specified frequencies
% # <matlab:web('doc_specs.html') Inserting the design specfications>
% # <matlab:web('doc_bounds.html') Computing bounds> over the Nichols chart that give constraint to the loop
% shaping of the open-loop $L(s) = P(s) K(s)$
% # Designing the feedback compensator $K(s)$ which satisfy the bounds
% # Designing the pre-filter (if needed) to shape the reference such that closed-loop
% sepcifications are satisfied



