% ==========================================
% --- RECEIVER: Audio over RF (Hybrid) ---
% ==========================================
clear; clc; close all;
set(0,'DefaultFigureVisible','on');

%% --- USER CONFIGURATION ---
mod_list = {'BPSK', 'QPSK', '8-QAM', '16-QAM'};
[mod_idx, tf] = listdlg('PromptString', 'Select a Modulation Scheme:', ...
                        'SelectionMode', 'single', ...
                        'ListString', mod_list, ...
                        'Name', 'TX/RX Config', ...
                        'ListSize', [200, 100]);
if tf == 0
    error('Configuration canceled. Script stopped.');
end
mod_scheme = mod_list{mod_idx};
disp(['Selected Modulation: ', mod_scheme]);
filter_type = 'RRC25'; 
% Choose: 'Rect', 'RRC25', 'RRC50', 'RRC75'

% Dynamically set Sync Bit length
if strcmp(mod_scheme, '8-QAM')
    N = 63; 
else
    N = 64; 
end
samplesPerSymbol = 1;

h_pn = comm.PNSequence('Polynomial', [6 5 0], 'InitialConditions', [0 0 0 0 0 1], ...
    'VariableSizeOutput', true, 'MaximumOutputSize', [N, 1]);
syncronization_bits = h_pn(N);

%% 1. Read LabVIEW Data and Filter
RX_ss_tmp = csvread('Acquired_data.csv',0,0); 
RX_ss = complex(RX_ss_tmp(:,1), RX_ss_tmp(:,2)); 

switch filter_type
    case 'Rect',  FIR_coeff = ones(samplesPerSymbol, 1) / samplesPerSymbol;
    case 'RRC25', FIR_coeff = rcosdesign(0.25, 8, samplesPerSymbol, 'sqrt').';
    case 'RRC50', FIR_coeff = rcosdesign(0.50, 8, samplesPerSymbol, 'sqrt').';
    case 'RRC75', FIR_coeff = rcosdesign(0.75, 8, samplesPerSymbol, 'sqrt').';
end

XX_signal = conv(FIR_coeff, RX_ss);
RX_signal = XX_signal / max([real(XX_signal); imag(XX_signal)]);

%% 2. EYE Diagram & Sampling Offset
% To see the "eye," we need to overlay multiple traces of 2 symbol periods
eye_len = samplesPerSymbol * 2; 
eye_frame_len = floor(length(RX_signal)/samplesPerSymbol) - 1;

I_eye = zeros(eye_len, eye_frame_len);
for i = 1:eye_frame_len
    % Slide by 1 symbol at a time to capture all transitions
    idx_start = (i-1)*samplesPerSymbol + 1;
    I_eye(:,i) = real(RX_signal(idx_start : idx_start + eye_len - 1));
end

% --- Plotting with Dark Theme ---
figure('Color', [0.1 0.1 0.1]);
% Plot the first 200 traces with some transparency (0.3) for a professional look
p = plot(I_eye(:, 1:min(200, eye_frame_len)), 'Color', [0.38 0.55 0.95, 0.3]); 
title(['I-Eye Diagram | ', mod_scheme], 'Color', 'w');
xlabel('Samples (2 Symbol Periods)', 'Color', 'w');
ylabel('Amplitude', 'Color', 'w');

% Style the axes
ax = gca;
ax.Color = [0.05 0.05 0.05];
ax.XColor = 'w'; ax.YColor = 'w';
ax.GridColor = [0.4 0.4 0.4];
grid on;
ylim([-1.2 1.2]); % Set limits to see the full BPSK/QAM range

% --- Sampling Offset Calculation ---
% We look for the point in the symbol period where the variance is highest 
% (the widest part of the eye)
eye_var = zeros(samplesPerSymbol, 1);
for i = 1:eye_frame_len
    for j = 1:samplesPerSymbol
        eye_var(j) = eye_var(j) + I_eye(j, i)^2;
    end
end
eye_var = eye_var / eye_frame_len;
[~, eye_offset] = max(eye_var);

disp(['Calculated Optimal Sampling Offset: ', num2str(eye_offset)]);

%% 3. Downsampling
Symbols = zeros(floor(length(RX_signal)/samplesPerSymbol),1);
for i=1:length(Symbols)-1
    Symbols(i,1) = RX_signal((i-1)*samplesPerSymbol+eye_offset,1);
end
avg_pwr = mean(abs(Symbols).^2);
Symbols = Symbols / sqrt(avg_pwr);

%% 4 & 5. PLL, Demodulation, and Synchronization
tx_ref_bits = csvread('Transmitted_ref.csv', 0, 0);

if strcmp(mod_scheme, 'BPSK')
    K1_PLL = 0.0313; K2_PLL = 2.49e-4;
    unwrap_phi_array = zeros(size(Symbols));
    phi_array = zeros(size(Symbols));
    filt_phi_array = zeros(size(Symbols));
    phi_array(1:2) = atan(imag(Symbols(1:2))./real(Symbols(1:2)));
    unwrap_phi_array(1:2) = phi_array(1:2);
    for i=3:length(Symbols)
        phi = atan(imag(Symbols(i))/real(Symbols(i)));
        old_phi = phi_array(i-1);
        if (phi-old_phi < -pi/2), freq = pi+phi-old_phi;
        elseif (phi-old_phi > pi/2), freq = -pi+phi-old_phi;
        else, freq = phi-old_phi;
        end
        unwrap_phi_array(i) = unwrap_phi_array(i-1)+freq;
        phi_array(i) = phi;
        filt_phi_array(i) = (2-K1_PLL-K2_PLL)*filt_phi_array(i-1) -(1-K1_PLL)*filt_phi_array(i-2) ...
            +(K1_PLL+K2_PLL)*unwrap_phi_array(i-1) -(K1_PLL)*unwrap_phi_array(i-2);
        Symbols(i) = Symbols(i)*complex(cos(filt_phi_array(i)),-sin(filt_phi_array(i)));
    end
    bits_rec_bipolar = sign(real(Symbols));
    corval = zeros(length(bits_rec_bipolar)-64,1);
    for i=1:length(corval)
        corval(i,1) = sum(syncronization_bits(1:64) .* bits_rec_bipolar(i:i+63));
    end
    [~, start_bit_ind] = max(abs(corval));
    ss = sign(corval(start_bit_ind));        
    bits_rec = (ss * bits_rec_bipolar + 1) / 2;
else
    switch mod_scheme
        case 'QPSK',   M = 4;  bps = 2;
        case '8-QAM',  M = 8;  bps = 3; 
        case '16-QAM', M = 16; bps = 4;
    end
    sync_pad_len = mod(-length(syncronization_bits), bps);
    sync_sym = qammod([syncronization_bits; zeros(sync_pad_len, 1)], M, 'InputType', 'bit', 'UnitAveragePower', true);
    
    K1_PLL = 0.05; K2_PLL = 0.001; 
    phase_est = 0; freq_est = 0;
    Symbols_sync = zeros(size(Symbols));
    for k = 1:length(Symbols)
        sym_rot = Symbols(k) * exp(-1i * phase_est);
        Symbols_sync(k) = sym_rot;
        dec_bit = qamdemod(sym_rot, M, 'UnitAveragePower', true);
        sym_ideal = qammod(dec_bit, M, 'UnitAveragePower', true);
        phase_err = angle(sym_rot * conj(sym_ideal));
        freq_est = freq_est + K2_PLL * phase_err;
        phase_est = phase_est + K1_PLL * phase_err + freq_est;
    end
    corval = zeros(length(Symbols_sync) - length(sync_sym), 1);
    for i = 1:length(corval)
        corval(i) = abs(sum(conj(sync_sym) .* Symbols_sync(i : i+length(sync_sym)-1)));
    end
    [~, start_sym_ind] = max(corval);
    preamble_rx = Symbols_sync(start_sym_ind : start_sym_ind + length(sync_sym) - 1);
    phase_ambiguity = angle(sum(preamble_rx .* conj(sync_sym)));
    Symbols_corrected = Symbols_sync * exp(-1i * phase_ambiguity);
    
    num_symbols_expected = length(tx_ref_bits) / bps;
    rx_packet_sym = Symbols_corrected(start_sym_ind : start_sym_ind + num_symbols_expected - 1);
    bits_rec = qamdemod(rx_packet_sym, M, 'OutputType', 'bit', 'UnitAveragePower', true);
    start_bit_ind = 1; 
end

%% 6. Extract Header and Audio Data
if strcmp(mod_scheme, 'BPSK')
    idx = start_bit_ind + N;     
else
    idx = N + 1;                 
end

% --- Extract Audio Length (24 bits majority voting) ---
len_rep = bits_rec(idx : idx + 72 - 1);
len_majority = zeros(24,1);
for i = 1:24
    bits_triplet = len_rep(3*i-2 : 3*i);
    len_majority(i) = sum(bits_triplet) >= 2;
end
rx_num_samples = bin2dec(char(len_majority + '0')');
idx = idx + 72;

% --- Extract Sample Rate (24 bits majority voting) ---
fs_rep = bits_rec(idx : idx + 72 - 1);
fs_majority = zeros(24,1);
for i = 1:24
    bits_triplet = fs_rep(3*i-2 : 3*i);
    fs_majority(i) = sum(bits_triplet) >= 2;
end
rx_fs = bin2dec(char(fs_majority + '0')');
idx = idx + 72;

disp(['Header Info: ', num2str(rx_num_samples), ' samples @ ', num2str(rx_fs), ' Hz']);

% --- Extract Audio Data ---
audio_bit_len = rx_num_samples * 6;
rx_audio_bits = bits_rec(idx : idx + audio_bit_len - 1);

%% 7. Bit Error Rate (BER) Calculation
if strcmp(mod_scheme, 'BPSK')
    total_packet_length = N + 72 + 72 + audio_bit_len; 
    received_packet = bits_rec(start_bit_ind : start_bit_ind + total_packet_length - 1);
else
    received_packet = bits_rec(1:length(tx_ref_bits));
end
Nof_err_bits = sum(abs(tx_ref_bits - received_packet)); 
BER = Nof_err_bits / length(tx_ref_bits);
disp(['BER: ', num2str(BER)]);

%% 8. Reconstruct and Play Audio
audio_bin_rec = reshape(rx_audio_bits, 6, []).';
audio_uint6 = uint8(bin2dec(char(audio_bin_rec + '0')));
audio_final = (double(audio_uint6) / 31.5) - 1.0;

% % Scaling: Convert 0:255 back to -1.0:1.0
% audio_final = (double(audio_uint8) / 127.5) - 1.0;

figure; plot(audio_final); title(['Recovered Audio | BER: ', num2str(BER)]);
if rx_fs > 0 && rx_fs <= 192000
    sound(audio_final, rx_fs);
else
    disp('Sample rate out of bounds, playback skipped.');
end

%% 9. Save Reconstructed Audio to File
output_filename = 'Received_Audio_Output.wav';

% Check if audio data exists before saving
if ~isempty(audio_final)
    % Ensure data is in the correct range for audioread (-1 to 1)
    % We already did this in Step 8, but we clip to be safe
    audio_to_save = max(min(audio_final, 1), -1);
    
    audiowrite(output_filename, audio_to_save, rx_fs);
    
    fprintf('--- Process Complete ---\n');
    fprintf('Audio saved as: %s\n', output_filename);
    fprintf('Path: %s\n', fullfile(pwd, output_filename));
else
    disp('Error: No audio data available to save.');
end

%% --- PERFORMANCE ANALYSIS ---
% Choose the correct synchronized symbols based on modulation
if strcmp(mod_scheme, 'BPSK')
    sym_for_analysis = Symbols(start_bit_ind : start_bit_ind + length(received_packet) - 1);
else
    sym_for_analysis = rx_packet_sym;
end

% 1. Signal-to-Noise Ratio (SNR) Estimation
% Noise power = Mean Square Error between received symbols and ideal constellation points
M_val = 2; % default for BPSK
if ~strcmp(mod_scheme, 'BPSK')
    switch mod_scheme
        case 'QPSK', M_val = 4;
        case '8-QAM', M_val = 8;
        case '16-QAM', M_val = 16;
    end
end

ideal_syms = qammod(qamdemod(sym_for_analysis, M_val, 'UnitAveragePower', true), M_val, 'UnitAveragePower', true);
noise = sym_for_analysis - ideal_syms;
sig_pwr = mean(abs(ideal_syms).^2);
noise_pwr = mean(abs(noise).^2);

estimated_SNR = 10 * log10(sig_pwr / noise_pwr);
bits_per_sym = log2(M_val);
estimated_EbNo = estimated_SNR - 10 * log10(bits_per_sym);

fprintf('\n--- Performance Metrics ---\n');
fprintf('Estimated SNR: %.2f dB\n', estimated_SNR);
fprintf('Estimated Eb/N0: %.2f dB\n', estimated_EbNo);
fprintf('Measured BER: %e\n', BER);

%% --- VISUALIZATION ---
% Constellation Diagram
figure('Color', 'w', 'Name', 'Constellation Analysis');
plot(sym_for_analysis, '.', 'Color', [0.8 0.2 0.2], 'MarkerSize', 6); hold on;
plot(ideal_syms, 'k+', 'LineWidth', 2, 'MarkerSize', 10);
grid on; axis square;
title(['Constellation: ', mod_scheme, ' (SNR: ', num2str(round(estimated_SNR,2)), ' dB)']);
xlabel('In-Phase'); ylabel('Quadrature');
legend('Received Symbols', 'Ideal Points');

% Power Spectral Density (PSD)
figure('Color', 'w', 'Name', 'Frequency Domain Analysis');
periodogram(RX_signal, rectwin(length(RX_signal)), length(RX_signal), 1, 'centered');
title(['Power Spectral Density of Received Signal (Normalized)']);