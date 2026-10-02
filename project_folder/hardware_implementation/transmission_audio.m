% ==========================================
% --- TRANSMITTER: Audio over RF (Hybrid) ---
% ==========================================
clear; clc;

%% --- USER CONFIGURATION ---
mod_list = {'BPSK', 'QPSK', '8-QAM', '16-QAM'};
[mod_idx, tf] = listdlg('PromptString', 'Select a Modulation Scheme:', ...
                        'SelectionMode', 'single', ...
                        'ListString', mod_list, ...
                        'Name', 'TX/RX Config', ...
                        'ListSize', [200, 100]);
if tf == 0, error('Configuration canceled.'); end
mod_scheme = mod_list{mod_idx};
filter_type = 'Rect'; 

% Dynamic Sync Bit length
if strcmp(mod_scheme, '8-QAM'), N = 63; else, N = 64; end
samplesPerSymbol = 1;

%% 1. Read and Process Audio
% Load audio file (Ensure the file exists in your path)
[audio_data, fs] = audioread('../files/input_audio.mpeg'); 

% Convert to Mono if Stereo
if size(audio_data, 2) > 1
    audio_data = mean(audio_data, 2); 
end

% --- Downsample Audio ---
target_fs = 3001; % choose 8000 / 16000 / 22050  (sound play on above 3000)

if fs ~= target_fs
    audio_data = resample(audio_data, target_fs, fs);
    fs = target_fs; % update for header transmission
end

% Downsample or truncate to keep packet size reasonable for RF simulation
% (Optional: Sending 1 second of audio)
max_duration = 0.5; % seconds
max_samples = fs * max_duration;
if length(audio_data) > max_samples
    audio_data = audio_data(1:max_samples);
end

% --- 6-bit Quantization ---
audio_6bit = uint8((audio_data + 1) * 31.5); % Range: 0–63

num_samples = length(audio_6bit);

% Convert to binary (6 bits instead of 8)
audio_bin = dec2bin(audio_6bit, 6).';
audio_bits = audio_bin(:) - '0';

figure; plot(audio_data); title(' Audio ');
%% 2. Create Header (Sync Bits + Audio Length + Sample Rate)
h_pn = comm.PNSequence('Polynomial', [6 5 0], 'InitialConditions', [0 0 0 0 0 1], ...
    'VariableSizeOutput', true, 'MaximumOutputSize', [N, 1]);
sync_bits = h_pn(N);

% We need to send the number of samples so the receiver knows when to stop
len_bits = dec2bin(num_samples, 24).'; % Support up to 16M samples
len_bits = len_bits(:) - '0';

% Optional: Send FS (Sample Rate) so the receiver plays it back correctly
fs_bits = dec2bin(fs, 24).'; 
fs_bits = fs_bits(:) - '0';

%% 3. Build the Packet
% Apply the same 3x repetition coding you used for image headers
len_bits_rep = reshape(repmat(len_bits.', 3, 1), [], 1);
fs_bits_rep = reshape(repmat(fs_bits.', 3, 1), [], 1);

tx_packet = [sync_bits; len_bits_rep; fs_bits_rep; audio_bits];

%% 4. Modulation & Pulse Shaping
switch filter_type
    case 'Rect',  FIR_coeff = ones(samplesPerSymbol, 1) / samplesPerSymbol;
    case 'RRC25', FIR_coeff = rcosdesign(0.25, 8, samplesPerSymbol, 'sqrt').';
    case 'RRC50', FIR_coeff = rcosdesign(0.50, 8, samplesPerSymbol, 'sqrt').';
    case 'RRC75', FIR_coeff = rcosdesign(0.75, 8, samplesPerSymbol, 'sqrt').';
end

if strcmp(mod_scheme, 'BPSK')
    bits_I = 2 * (tx_packet - 0.5); 
    sig_ipt_I = upsample(bits_I, samplesPerSymbol);
    tx_signal = conv(FIR_coeff, sig_ipt_I);
    tx_signal = tx_signal / max(abs(tx_signal));
    
    csvwrite('Transmitted_ref.csv', tx_packet); 
    csvwrite('Transmitted_data.csv', tx_signal);
else
    switch mod_scheme
        case 'QPSK',   M = 4;  bps = 2;
        case '8-QAM',  M = 8;  bps = 3; 
        case '16-QAM', M = 16; bps = 4;
    end
    pad_len = mod(-length(tx_packet), bps);
    tx_packet_padded = [tx_packet; zeros(pad_len, 1)];
    tx_sym = qammod(tx_packet_padded, M, 'InputType', 'bit', 'UnitAveragePower', true);
    sig_upsampled = upsample(tx_sym, samplesPerSymbol);
    
    tx_signal = conv(FIR_coeff, sig_upsampled);
    tx_signal = tx_signal / max(max(abs(real(tx_signal))), max(abs(imag(tx_signal))));
    
    csvwrite('Transmitted_ref.csv', tx_packet_padded); 
    csvwrite('Transmitted_data.csv', [real(tx_signal), imag(tx_signal)]); 
end

disp('--- Audio TX Complete ---');
disp(['Samples Sent: ', num2str(num_samples), ' at ', num2str(fs), ' Hz']);