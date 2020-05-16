function [res] = nicresp(G,frequency)
%NICRESP Computes frequenct response in Nichols form for general object


w = frequency;

switch class(G)
    %case 'qfr'
    %    res = G.nicresp(w);
    case {'tf','zpk','ss','frd'}
        Bfr = squeeze(freqresp(G,w)).';
        res = c2n(Bfr.',-180);
    %case 'qctrl'
    %    res = nfr(G,w);
    case 'double'
        %if ~isscalar(B), error('a numeric value must be a scalar'); end
        res = repmat(c2n(G),length(w),1);
    otherwise
        error('unsupported input class')
end



end

