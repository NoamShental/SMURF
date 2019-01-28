function [bact_freq, keep_col] = ml_em_iterative(A_L2, y_vec, Config, bact_freq)

keep_col = 1:size(A_L2,2);

if ~exist('bact_freq','var') || isempty(bact_freq)
    bact_freq = A_L2'*y_vec;
    bact_freq = bact_freq/sum(bact_freq);
end    


tStart = tic; tic
for jj = 1:Config.numIter
    
    % Estimate theta
    theta_i = A_L2*bact_freq;
    
    % Reweight the counts
    r_weighted = y_vec./(theta_i+eps);
    
    % Update the the factor
    bact_factor = r_weighted'*A_L2;

    
    % Calculate the error
    L1_error = abs(1-bact_factor)*bact_freq;

    % Update the frequency
    bact_freq = bact_freq.*bact_factor';
    
    % Check if rached tolerance
    if L1_error < Config.tol
        break
    end
    
    tTmp = toc;
    if tTmp > 60
        % Remove bacterias with zero frequency
        remove_ind = find(bact_freq < 1e-10);
        A_L2(:,remove_ind) = [];
        keep_col(remove_ind) = [];
        bact_freq(remove_ind) = [];
        
        % Print
        if Config.verbose
            disp(['Iter:' num2str(jj) '. Error reduction of X (L1 norm): ' num2str(L1_error)])
        end
        tic
    end
    
end

if Config.verbose
    disp(['Total iterations time: ' num2str(toc(tStart))])
end

% Hard thresholding
bact_freq(bact_freq<1e-10) = 0;
bact_freq = bact_freq/sum(bact_freq);

