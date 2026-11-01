% Richard Harris
% 10/10/2026
% rharr124@asu.edu
% EEE 596 Adaptive Signal Processing
% RLS Adaptive Equalizer Function

function [y] = rls_eq(x, d, L, lambda, delta)
    N = numel(x);
    y = complex(zeros(N, 1));
    bhat = complex(zeros(L, 1));
    k = complex(zeros(L, 1)); % gain vector
    P = (1/delta)*complex(eye(L, L)); % inverse of the autocorrelation matrix
    
    % RLS adaptive loop
    mse_est = zeros(N, 1);
    for i=1:N
        % compute the output of the adaptive system
        % y[n] = b0x[n] + b1x[n-1] + b2x[n-2] + ...
        if i>=L
            x_block = x(i:-1:i-L+1); % most recent sample is index i, aligns with b0
        else
            x_block = [x(i:-1:1); complex(zeros(L-i,1))]; % assume 0s before first sample
        end
        yhat = bhat.' * x_block; % .' is standard transpose (not Hermitian)
        y(i) = yhat;
    
        % compute the error and and mean squared error estimate
        e = d(i) - yhat;
        mse_est(i) = abs(e)*abs(e);

        % calculate the gain vector
        k = ((1/lambda)*P*x_block) ./ (1 + (1/lambda)*x_block.'*P*x_block);
    
        % update the coefficients
        % in vector form: b(n+1) = b(n) + k(n)*e
        bhat = bhat + k*e;

        % update the inverse of the autocorrelation matrix
        P = (1/lambda)*P - (1/lambda)*k*x_block.'*P;
    end
    
    % Plot MSE in dB vs iteration number
    % mse_est_dB = 10*log10(mse_est);
    % plot(mse_est_dB);
    % xlabel("Iteration");
    % ylabel("MSE (dB)");
    % title(sprintf("MSE (dB) vs Iteration number, L=%d, mu=%.6f", L, mu));
end
