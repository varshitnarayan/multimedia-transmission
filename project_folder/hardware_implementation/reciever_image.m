%% --- DEBUGGING RECEIVER: BPSK Image over RF ---
clear; clc; close all;
set(0,'DefaultFigureVisible','on');
%% 1. Parameters & Data Loading
N = 64; 
samplesPerSymbol = 2; % Matches your TX constraint
mod_scheme = 'BPSK';
rrc_alpha = 0.0; % Set to 0.0 for Rect filter, >0 for RRC
try
    RX_ss_tmp = csvread('Acquired_data_rrc75.csv',0,0);
    tx_ref_bits = csvread('Transmitted_ref.csv', 0, 0);
    RX_ss = complex(RX_ss_tmp(:,1), RX_ss_tmp(:,2));
    fprintf('Success: Data files loaded.\n');
catch
    error('File Error: Ensure Acquired_data.csv and Transmitted_ref.csv exist.');
end

%% 2. Filtering & Normalization (ADAPTIVE LOGIC)
if rrc_alpha == 0
    % Time-domain Rectangular Matched Filter (Moving Average)
    FIR_coeff = ones(samplesPerSymbol, 1) / samplesPerSymbol;
    fprintf('Using Rectangular (Boxcar) Matched Filter.\n');
else
    % Standard Root Raised Cosine Filter
    FIR_coeff = rcosdesign(rrc_alpha, 8, samplesPerSymbol, 'sqrt').';
    fprintf('Using RRC Matched Filter with Alpha = %.2f.\n', rrc_alpha);
end

XX_signal = conv(FIR_coeff, RX_ss);
RX_signal = XX_signal / max(abs(XX_signal));

% --- DEBUG PLOT: Time Domain ---
figure('Name', 'Debug: Raw vs Filtered');
subplot(2,1,1); plot(real(RX_ss(1:500))); title('Raw Received Signal (I)');
subplot(2,1,2); plot(real(RX_signal(1:500))); title('Filtered Signal (I)');

%% 3. Eye Diagram & Sampling
eye_len = samplesPerSymbol;
eye_frame_len = floor(length(RX_signal)/eye_len);
I_eye = reshape(real(RX_signal(1:eye_frame_len*eye_len)), eye_len, eye_frame_len);
figure('Name', 'Debug: Eye Diagram');
plot(I_eye(:, 1:min(200, end))); title('Eye Diagram (I-channel)');
grid on;
% Best sampling instant
eye_var = mean(I_eye.^2, 2);
[~, eye_offset] = max(eye_var);
fprintf('Debug: Best sampling offset found at index %d\n', eye_offset);
% Downsampling
Symbols = RX_signal(eye_offset:samplesPerSymbol:end);

%% 4. PLL (Phase/Frequency Recovery)
K1_PLL = 0.03; 
K2_PLL = 0.0035;
unwrap_phi_array = zeros(size(Symbols));
phi_array = zeros(size(Symbols));
filt_phi_array = zeros(size(Symbols));
% Initial phase
phi_array(1:2) = atan2(imag(Symbols(1:2)), real(Symbols(1:2)));
unwrap_phi_array(1:2) = phi_array(1:2);
for i=3:length(Symbols)
    phi = atan2(imag(Symbols(i)), real(Symbols(i)));
    old_phi = phi_array(i-1);
    
    delta_phi = phi - old_phi;
    if delta_phi < -pi/2,     freq = pi + delta_phi;
    elseif delta_phi > pi/2,  freq = -pi + delta_phi;
    else,                     freq = delta_phi;
    end
    
    unwrap_phi_array(i) = unwrap_phi_array(i-1) + freq;
    phi_array(i) = phi;
    
    filt_phi_array(i) = (2-K1_PLL-K2_PLL)*filt_phi_array(i-1) ...
        -(1-K1_PLL)*filt_phi_array(i-2) ...
        +(K1_PLL+K2_PLL)*unwrap_phi_array(i-1) ...
        - K1_PLL*unwrap_phi_array(i-2);
            
    Symbols(i) = Symbols(i) * exp(-1j * filt_phi_array(i));
end
figure('Name', 'Debug: PLL Constellation');
scatter(real(Symbols), imag(Symbols), '.'); 
title('Constellation After PLL'); grid on; xlim([-1.5 1.5]); ylim([-1.5 1.5]);

%% 5. Synchronization & Correlation
h_pn = comm.PNSequence('Polynomial', [6 5 0], 'InitialConditions', [0 0 0 0 0 1], ...
    'VariableSizeOutput', true, 'MaximumOutputSize', [N, 1]);
sync_bits = h_pn(N);
sync_bipolar = 2*sync_bits - 1;
bits_rec_bipolar = sign(real(Symbols));
corval = zeros(length(bits_rec_bipolar)-N, 1);
for i=1:length(corval)
    corval(i) = sum(sync_bipolar .* bits_rec_bipolar(i:i+N-1));
end
figure('Name', 'Debug: Correlation');
plot(corval); title('Synchronization Correlation');
[peak, start_bit_ind] = max(abs(corval));
fprintf('Debug: Max Correlation Peak = %f at index %d\n', peak, start_bit_ind);
if peak < (N * 0.6)
    warning('Sync Weak! Peak is only %d/%d. Header may be corrupt.', peak, N);
end
ss = sign(corval(start_bit_ind));        
bits_rec = (ss * bits_rec_bipolar + 1) / 2;

%% 6. Header & Image Reconstruction
try
    idx = start_bit_ind + N; 
    rows_rep = bits_rec(idx : idx + 47);
    rows_bits = zeros(16,1);
    for i = 1:16
        rows_bits(i) = sum(rows_rep(3*i-2 : 3*i)) >= 2;
    end
    rx_rows = bin2dec(char(rows_bits + '0')');
    idx = idx + 48;
    cols_rep = bits_rec(idx : idx + 47);
    cols_bits = zeros(16,1);
    for i = 1:16
        cols_bits(i) = sum(cols_rep(3*i-2 : 3*i)) >= 2;
    end
    rx_cols = bin2dec(char(cols_bits + '0')');
    idx = idx + 48;
    fprintf('Debug: Decoded Header -> %d x %d\n', rx_rows, rx_cols);
    img_bit_len = rx_rows * rx_cols * 8;
    rx_img_bits = bits_rec(idx : idx + img_bit_len - 1);
    
    img_bin_matrix = reshape(rx_img_bits, 8, []).'; 
    img_dec = bin2dec(char(img_bin_matrix + '0'));
    final_image = uint8(reshape(img_dec, rx_rows, rx_cols));
    figure('Name', 'Final Result');
    imshow(final_image);
    title(sprintf('Received Image %dx%d', rx_rows, rx_cols));
catch ME
    fprintf('CRITICAL ERROR during Header/Image parsing: %s\n', ME.message);
end

%% 7. BER Calculation
total_expected = length(tx_ref_bits);
if length(bits_rec) >= (start_bit_ind + total_expected - 1)
    received_packet = bits_rec(start_bit_ind : start_bit_ind + total_expected - 1);
    errors = sum(abs(tx_ref_bits - received_packet));
    fprintf('--- STATISTICS ---\n');
    fprintf('BER: %f\n', errors / total_expected);
    fprintf('Bit Errors: %d / %d\n', errors, total_expected);
else
    fprintf('Could not calculate BER: Received stream shorter than reference.\n');
end 

%% --- MANUAL ALIGNMENT CHECK ---
try
    ref_sig = csvread('Transmitted_data.csv');
    if size(ref_sig,2) > 1, ref_sig = ref_sig(:,1); end
    rx_raw = csvread('Acquired_data.csv'); 
    rx_sig = complex(rx_raw(:,1), rx_raw(:,2));
    rx_real = real(rx_sig);
    [xc, lags] = xcorr(rx_real, ref_sig);
    [~, max_idx] = max(abs(xc));
    best_lag = lags(max_idx);
    fprintf('The transmitted signal starts at index: %d\n', best_lag);
    figure('Name', 'Manual Signal Match Check');
    subplot(2,1,1); plot(ref_sig(1:min(2000, end))); title('Original Transmitted Signal'); grid on;
    subplot(2,1,2);
    if best_lag > 0
        matched_portion = rx_real(best_lag : min(best_lag + length(ref_sig)-1, length(rx_real)));
        plot(matched_portion(1:min(2000, end)));
        title(['Acquired Data (Starting at index ', num2str(best_lag), ')']);
    end
    grid on;
catch
    fprintf('Manual alignment check skipped (Files missing).\n');
end

%% 8. Extended Performance Metrics
if exist('start_bit_ind', 'var') && start_bit_ind > 0
    eval_syms = Symbols(start_bit_ind : min(start_bit_ind + length(bits_rec_bipolar) - 1, end));
    ideal_pts = ss * sign(real(eval_syms)); 
    errors_vec = eval_syms - ideal_pts;
    evm_rms = sqrt(mean(abs(errors_vec).^2)) / mean(abs(ideal_pts)) * 100;
    sig_pwr = mean(abs(ideal_pts).^2);
    noise_pwr = mean(abs(errors_vec).^2);
    snr_est = 10 * log10(sig_pwr / noise_pwr);
    
    fprintf('\n--- RF SIGNAL QUALITY ---\n');
    fprintf('RMS EVM: %.2f%%\n', evm_rms);
    fprintf('Estimated SNR: %.2f dB\n', snr_est);
    
    figure('Name', 'Metric: EVM Analysis');
    hold on;
    scatter(real(eval_syms), imag(eval_syms), '.', 'MarkerEdgeAlpha', 0.3);
    plot([-1 1], [0 0], 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    title(['EVM: ', num2str(evm_rms, '%.2f'), '% | Est. SNR: ', num2str(snr_est, '%.2f'), ' dB']);
    grid on; hold off;
end