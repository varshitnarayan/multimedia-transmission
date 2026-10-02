%% --- AUDIO RECEIVER: BPSK matched to 6-bit TX ---
clear; clc; close all;
%% 1. Parameters (MUST MATCH TRANSMITTER)
N = 64; 
samplesPerSymbol = 1; % Changed to 1 to match your TX
mod_scheme = 'BPSK';
alpha = 0.50; % Set to 0.0 for Rect filter, >0 for RRC

% Load data
try
    RX_ss_tmp = csvread(['acquire_50_audio' ...
        '.csv'], 0, 0);
    tx_ref_bits = csvread('Transmitted_ref_aud_rrc50.csv', 0, 0);
    RX_ss = complex(RX_ss_tmp(:,1), RX_ss_tmp(:,2));
    fprintf('Success: Data files loaded.\n');
catch
    error('File Error: Ensure Acquired_data.csv and Transmitted_ref.csv exist.');
end

%% 2. Filtering & Recovery (ADAPTIVE LOGIC)
if alpha == 0
    % Time-domain Rectangular Matched Filter (Moving Average)
    FIR_coeff = ones(samplesPerSymbol, 1) / samplesPerSymbol;
    fprintf('Using Rectangular (Boxcar) Matched Filter.\n');
else
    % Standard Root Raised Cosine Filter
    FIR_coeff = rcosdesign(alpha, 8, samplesPerSymbol, 'sqrt').';
    fprintf('Using RRC Matched Filter with Alpha = %.2f.\n', alpha);
end

RX_signal = conv(FIR_coeff, RX_ss);
RX_signal = RX_signal / max(abs(RX_signal));

% Best sampling instant
eye_len = samplesPerSymbol;
if eye_len > 1
    eye_frame_len = floor(length(RX_signal)/eye_len);
    I_eye = reshape(real(RX_signal(1:eye_frame_len*eye_len)), eye_len, eye_frame_len);
    [~, eye_offset] = max(mean(I_eye.^2, 2));
    Symbols = RX_signal(eye_offset:samplesPerSymbol:end);
else
    Symbols = RX_signal; % If SPS=1, no downsampling needed
end

%% 3. PLL (Using your working "aggressive" constants)
K1_PLL = 0.03; K2_PLL = 0.0035;
unwrap_phi_array = zeros(size(Symbols));
phi_array = zeros(size(Symbols));
filt_phi_array = zeros(size(Symbols));
phi_array(1:2) = atan2(imag(Symbols(1:2)), real(Symbols(1:2)));
unwrap_phi_array(1:2) = phi_array(1:2);
for i=3:length(Symbols)
    phi = atan2(imag(Symbols(i)), real(Symbols(i)));
    delta_phi = phi - phi_array(i-1);
    if delta_phi < -pi/2, freq = pi + delta_phi;
    elseif delta_phi > pi/2, freq = -pi + delta_phi;
    else, freq = delta_phi;
    end
    unwrap_phi_array(i) = unwrap_phi_array(i-1) + freq;
    phi_array(i) = phi;
    filt_phi_array(i) = (2-K1_PLL-K2_PLL)*filt_phi_array(i-1) - (1-K1_PLL)*filt_phi_array(i-2) + (K1_PLL+K2_PLL)*unwrap_phi_array(i-1) - K1_PLL*unwrap_phi_array(i-2);
    Symbols(i) = Symbols(i) * exp(-1j * filt_phi_array(i));
end

%% 4. Synchronization
h_pn = comm.PNSequence('Polynomial', [6 5 0], 'InitialConditions', [0 0 0 0 0 1], 'VariableSizeOutput', true, 'MaximumOutputSize', [N, 1]);
sync_bits = h_pn(N);
sync_bipolar = 2*sync_bits - 1;
bits_rec_bipolar = sign(real(Symbols));
corval = zeros(length(bits_rec_bipolar)-N, 1);
for i=1:length(corval)
    corval(i) = sum(sync_bipolar .* bits_rec_bipolar(i:i+N-1));
end
[~, start_bit_ind] = max(abs(corval));
ss = sign(corval(start_bit_ind));        
bits_rec = (ss * bits_rec_bipolar + 1) / 2;

%% 5. Header Extraction (Majority Voting)
idx = start_bit_ind + N;
% --- Decode num_samples (24 bits x 3 = 72 bits) ---
num_samples_rep = bits_rec(idx : idx + 71);
num_samples_bits = zeros(24,1);
for i = 1:24
    num_samples_bits(i) = sum(num_samples_rep(3*i-2 : 3*i)) >= 2;
end
rx_num_samples = bin2dec(char(num_samples_bits + '0')');
idx = idx + 72;
% --- Decode Sample Rate (24 bits x 3 = 72 bits) ---
fs_rep = bits_rec(idx : idx + 71);
fs_bits = zeros(24,1);
for i = 1:24
    fs_bits(i) = sum(fs_rep(3*i-2 : 3*i)) >= 2;
end
rx_fs = bin2dec(char(fs_bits + '0')');
idx = idx + 72;
fprintf('Audio Header: %d samples at %d Hz\n', rx_num_samples, rx_fs);

%% 6. Audio Data Reconstruction (6-bit)
total_audio_bits = rx_num_samples * 6;
rx_audio_bits = bits_rec(idx : idx + total_audio_bits - 1);
% Reshape into 6-bit chunks
audio_bin_matrix = reshape(rx_audio_bits, 6, []).'; 
audio_uint6 = bin2dec(char(audio_bin_matrix + '0'));
% Reverse the 6-bit Quantization: audio_6bit = uint8((audio_data + 1) * 31.5)
audio_out = (double(audio_uint6) / 31.5) - 1;

%% 7. Output
figure; plot(audio_out); title('Received Audio Waveform');
if rx_fs > 0 && rx_fs < 100000 % Safety check for sample rate
    sound(audio_out, rx_fs);
    audiowrite('received_audio.wav', audio_out, rx_fs);
end
% BER
total_expected = length(tx_ref_bits);
received_packet = bits_rec(start_bit_ind : min(start_bit_ind + total_expected - 1, end));
fprintf('BER: %f\n', sum(abs(tx_ref_bits(1:length(received_packet)) - received_packet)) / total_expected);

%% 8. Extended Performance Metrics
% EVM Calculation
ideal_symbols = ss * (2*bits_rec(start_bit_ind:start_bit_ind+length(Symbols)-start_bit_ind) - 1);
% Ensure vectors match in length for comparison
len_eval = min(length(Symbols(start_bit_ind:end)), length(ideal_symbols));
actual_syms = Symbols(start_bit_ind:start_bit_ind+len_eval-1);
ref_syms = ideal_symbols(1:len_eval);
evm_vec = abs(actual_syms - ref_syms).^2;
evm_rms = sqrt(mean(evm_vec)) / mean(abs(ref_syms)) * 100;
% Estimated SNR from EVM
snr_est = 20 * log10(1 / (evm_rms / 100));
fprintf('--- ADVANCED METRICS ---\n');
fprintf('RMS EVM: %.2f%%\n', evm_rms);
fprintf('Estimated SNR: %.2f dB\n', snr_est);
fprintf('------------------------\n');