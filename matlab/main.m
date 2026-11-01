% Richard Harris
% 10/10/2026
% rharr124@asu.edu
% EEE 596 Adaptive Signal Processing
% Project - LMS and RLS Equalizers in Optical Communications

% TODO: adaptive/normalized learning rate for LMS, rls parameters and
% why does it have larger BER for 0 and 2 SNR? mu=0.0005 for LMS achieves
% similar performance to RLS

close all; clear; clc;

rng(43);

% parameters
h = [0.74; -0.514; 0.37; 0.216; 0.062]; % channel impulse response
SNR_sweep = 0:2:16; % sweep SNR from 0 to 16 dB
N = 200000; % number of symbols
M = 2; % BPSK is modulation order 2
algs = ["none", "lms", "rls"];
L = 10; % adaptive filter order
mu = 0.00001; % LMS learning rate
lambda = 0.999; % RLS forgetting factor
delta = 1; % RLS autocorrelation initialization parameter, make larger for lower SNR

% BER loop for each SNR value
ber_sweep = zeros(numel(algs), numel(SNR_sweep));
for i = 1:numel(algs)
    alg = algs(i);
    for j = 1:numel(SNR_sweep)
        SNR = SNR_sweep(j);
    
        % transmitter
        data = randi([0,1], N, 1); % randomly generate BPSK symbol indices
        x = pskmod(data, M , 0); % BPSK modulated data
        
        % channel
        s = filter(h, 1, x); % FIR channel
        u = awgn(s, SNR); % add AWGN
        
        % receiver
        if alg=="none"
            y = u; % no equalizer
        elseif alg=="lms"
            y = lms_eq(u, x, L, mu);
        elseif alg=="rls"
            y = rls_eq(u, x, L, lambda, delta);
        else
            error("Algorithm not supported");
        end
        
        % received symbols
        r = pskdemod(y, M, 0); % demodulate to the closest BPSK constellation point
        bit_errors = sum(r ~= data);
        ber = bit_errors / N;
        ber_sweep(i,j) = ber;
        fprintf('Algorithm: %s, SNR = %d dB, BER = %.4f\n', alg, SNR, ber);
    end
    fprintf('\n');
end

% plot results
semilogy(SNR_sweep,ber_sweep(1,:),'-or','MarkerSize',6,'MarkerEdgeColor','red','MarkerFaceColor','red','DisplayName','Without Equalizer');
hold on;
semilogy(SNR_sweep,ber_sweep(2,:),'-ob','MarkerSize',6,'MarkerEdgeColor','blue','MarkerFaceColor','blue','DisplayName','LMS Equalizer');
semilogy(SNR_sweep,ber_sweep(3,:),'-og','MarkerSize',6,'MarkerEdgeColor','green','MarkerFaceColor','green','DisplayName','RLS Equalizer');
hold off;
ylim([1e-3,1]);
title('BER Performance of Equalizers in Optical Communication (Reproduced)');
xlabel('SNR (dB)');
ylabel('BER', "Rotation", 0);
grid on;
legend;
