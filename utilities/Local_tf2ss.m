function [a,b,c,d] = Local_tf2ss(num,den)
%LOCAL_TF2SS Converts two polynomial to state-space model (A,B,C,D)
%           Local_tf2ss(num,den) converts tf data to state-space model
%           without the control system toolbox, used for time domain
%           simulations. This is purely algebraic, sample time (if it
%           exists) is not retained.

    % Cast to enforce single precision rules
    if isa(num,'single') || isa(den,'single')
        nums = single(num);
        dens = single(den);
    else
        nums = num;
        dens = den;
    end
    %check if null system  (both numerator and denominator are empty)
    if isempty(nums) && isempty(dens)
        a = zeros(0,'like',nums);
        b = zeros(0,'like',nums);
        c = zeros(0,'like',nums);
        d = zeros(0,'like',nums);
    else
        assert(ismatrix(nums) && ismatrix(dens),'Inputs should be two-dimensional.');
        if(min(size(dens)) > 1)
            error('Denominator must be a row vector.');
        end
        denRow = dens(:).';
        % Index of first non zero element of denominator
        startIndexDen = find(denRow,1);
        % Denominator should not be zero or empty
        if isempty(startIndexDen)
            error('Denominator cannot be zero.');
        end
        % Strip denominator of leading zeros
        denStrip = denRow(startIndexDen(1):end);
        [mnum,nnum] = size(nums);
        nden = size(denStrip,2);
        % Check for proper numerator
        if (nnum > nden)
            if any(nums(:,1:(nnum - nden)) ~= 0,'all')
                error('Transfer function must be proper!');
            end
            % Try to strip leading zeros to make proper
            numStrip = nums(:,(nnum-nden+1):nnum);
        else
            % Pad numerator with leading zeroes, to make it have same number of
            % Columns as the denominator
            numStrip = [zeros(mnum,nden-nnum) nums];
        end

        % Normalize numerator and denominator such that first element of
        % Denominator is one
        numNormalized = numStrip./denStrip(1);
        denNormalized = denStrip./denStrip(1);
        if mnum == 0
            d = zeros(0,'like',numNormalized);
            c = zeros(0,'like',numNormalized);
        else
            d = numNormalized(:,1);
            c = numNormalized(:,2:nden) - numNormalized(:,1) * denNormalized(2:nden);
        end

        if nden == 1
            a = zeros(0,'like',numNormalized);
            b = zeros(0,'like',numNormalized);
            c = zeros(0,'like',numNormalized);
        else
            a = [-denNormalized(2:nden);eye(nden-2,nden-1)];
            b = eye(nden-1,1,'like',numNormalized);
        end
    end
end
