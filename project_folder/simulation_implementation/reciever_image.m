% ==========================================
% --- RECEIVER: Image over RF (Hybrid) ---
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

% Check if the user clicked "Cancel"
if tf == 0
    error('Configuration canceled. Script stopped.');
end

mod_scheme = mod_list{mod_idx};
disp(['Selected Modulation: ', mod_scheme]);

% mod_scheme = '8-QAM'; % Choose: 'BPSK', 'QPSK', '8-QAM', '16-QAM'
filter_type = 'RRC75'; % Choose: 'Rect', 'RRC25', 'RRC50', 'RRC75'

% Dynamically set Sync Bit length based on modulation
if strcmp(mod_scheme, '8-QAM')
    N = 63; % 63 bits perfectly divides by 3 bps (21 symbols)
else
    N = 64; % 64 bits perfectly divides by 1, 2, and 4 bps
end

samplesPerSymbol = 8;

h_pn = comm.PNSequence('Polynomial', [6 5 0], 'InitialConditions', [0 0 0 0 0 1], 'VariableSizeOutput', true, 'MaximumOutputSize', [N, 1]);
syncronization_bits = h_pn(N);

%% 1. Read LabVIEW Data and Filter
RX_ss_tmp = csvread('Acquired_data.csv',0,0); 
% USRP always provides I & Q on receive
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
% 1. Use 2 symbols per trace to see the full "eye"
eye_len = 2 * samplesPerSymbol; 
num_traces = floor(length(RX_signal) / samplesPerSymbol) - 1;

% Pre-allocate and use vectorization to fill traces
I_eye = zeros(eye_len, num_traces);
for i = 1:num_traces
    idx = (i-1)*samplesPerSymbol + 1;
    I_eye(:, i) = real(RX_signal(idx : idx + eye_len - 1));
end

% 2. Automated Offset Detection (Finding the widest part of the eye)
% We calculate variance across all traces; the peak is the optimal sampling point
eye_var_clean = var(I_eye, 0, 2); 
[~, eye_offset] = max(eye_var_clean(1:samplesPerSymbol));

% 3. Visualization
figure('Name', 'Signal Analysis');
subplot(2,1,1);
plot(I_eye(:, 1:min(200, end)), 'b'); % Plot first 200 traces
grid on; hold on;
xline(eye_offset, '--r', 'Optimal Sampling', 'LabelVerticalAlignment', 'bottom');
title(['I-Eye Diagram for ', mod_scheme]);

%% --- Revised Eye Diagram, Sampling & Constellation ---

% 1. Setup Eye Diagram Parameters
eye_len = 2 * samplesPerSymbol; % Display two symbols to see the full "eye"
num_traces = floor(length(RX_signal) / samplesPerSymbol) - 1;

% Pre-allocate and fill traces using the real part (I-channel)
I_eye = zeros(eye_len, num_traces);
for i = 1:num_traces
    idx = (i-1)*samplesPerSymbol + 1;
    I_eye(:, i) = real(RX_signal(idx : idx + eye_len - 1));
end

% 2. Automated Offset Detection
% We find where the variance is highest (widest eye opening)
eye_var_clean = var(I_eye, 0, 2); 
[~, eye_offset] = max(eye_var_clean(1:samplesPerSymbol));

% 3. Visualization
figure('Name', 'Communication System Analysis', 'Color', 'w');

% Top Subplot: Eye Diagram
subplot(2,1,1);
plot(I_eye(:, 1:min(300, end)), 'Color', [0.7 0.7 0.7]); % Light gray for background traces
hold on;
plot(I_eye(:, 1), 'b', 'LineWidth', 1.5); % Highlight one trace in blue
xline(eye_offset, '--r', 'Optimal Sampling', 'LabelVerticalAlignment', 'bottom');
grid on;
title(['I-Eye Diagram: ', mod_scheme]);
xlabel('Samples'); ylabel('Amplitude');

% Bottom Subplot: Constellation Diagram
subplot(2,1,2);
% Decimate the signal using the offset to get the actual symbol points
RX_symbols = RX_signal(eye_offset : samplesPerSymbol : end);

% Plotting the dots in RED ('r.')
plot(real(RX_symbols), imag(RX_symbols), 'r.', 'MarkerSize', 12); 

grid on; 
axis square;
title(['Received ', mod_scheme, ' Constellation (Red Dots)']);
xlabel('In-Phase (I)'); 
ylabel('Quadrature (Q)');

% Center the axis for a professional look
limit = max(abs([real(RX_symbols); imag(RX_symbols)])) * 1.5;
if limit > 0
    axis([-limit limit -limit limit]);
end

hold off;

%% 3. Downsampling
Symbols = zeros(floor(length(RX_signal)/samplesPerSymbol),1);
for i=1:length(Symbols)-1
    Symbols(i,1) = RX_signal((i-1)*samplesPerSymbol+eye_offset,1);
end

% --- THE FIX: Restore Average Power to 1 ---
% Undo the peak-normalization so 8-QAM decision boundaries align perfectly
avg_pwr = mean(abs(Symbols).^2);
Symbols = Symbols / sqrt(avg_pwr);
%% 4 & 5. PLL, Demodulation, and Synchronization
tx_ref_bits = csvread('Transmitted_ref.csv', 0, 0);

if strcmp(mod_scheme, 'BPSK')
    % ==========================================
    % --- ORIGINAL BPSK LOGIC (Untouched) ---
    % ==========================================
    K1_PLL = 0.0313;
    K2_PLL = 2.49e-4;
    unwrap_phi_array = zeros(size(Symbols));
    phi_array = zeros(size(Symbols));
    filt_phi_array = zeros(size(Symbols));

    phi_array(1:2) = atan(imag(Symbols(1:2))./real(Symbols(1:2)));
    unwrap_phi_array(1:2) = phi_array(1:2);

    for i=3:length(Symbols)
        phi = atan(imag(Symbols(i))/real(Symbols(i)));
        old_phi = phi_array(i-1);
        if (phi-old_phi < -pi/2)
            freq = pi+phi-old_phi;
        elseif (phi-old_phi > pi/2)
            freq = -pi+phi-old_phi;
        else
            freq = phi-old_phi;
        end
        unwrap_phi_array(i) = unwrap_phi_array(i-1)+freq;
        phi_array(i) = phi;
        
        filt_phi_array(i) = (2-K1_PLL-K2_PLL)*filt_phi_array(i-1) ...
            -(1-K1_PLL)*filt_phi_array(i-2) ...
            +(K1_PLL+K2_PLL)*unwrap_phi_array(i-1) ...
            -(K1_PLL)*unwrap_phi_array(i-2);
            
        Symbols(i) = Symbols(i)*complex(cos(filt_phi_array(i)),-sin(filt_phi_array(i)));
    end

    bits_rec_bipolar = sign(real(Symbols));
    corval = zeros(length(bits_rec_bipolar),1);

    for i=1:(length(bits_rec_bipolar)-64)
        corval(i,1) = sum(syncronization_bits(1:64) .* bits_rec_bipolar(i:i+63));
    end

    figure; plot(corval); title('Sync Correlation Peak (BPSK)');
    [peak, start_bit_ind] = max(abs(corval));

    ss = sign(corval(start_bit_ind));        
    bits_rec = (ss * bits_rec_bipolar + 1) / 2;
    
else
    % ==========================================
    % --- UPGRADED QAM LOGIC ---
    % ==========================================
    switch mod_scheme
        case 'QPSK',   M = 4;  bps = 2;
        case '8-QAM',  M = 8;  bps = 3; 
        case '16-QAM', M = 16; bps = 4;
    end
    
    sync_pad_len = mod(-length(syncronization_bits), bps);
    sync_sym = qammod([syncronization_bits; zeros(sync_pad_len, 1)], M, 'InputType', 'bit', 'UnitAveragePower', true);

    % Decision-Directed PLL
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

    % Complex Correlation
    corval = zeros(length(Symbols_sync) - length(sync_sym), 1);
    for i = 1:length(corval)
        corval(i) = abs(sum(conj(sync_sym) .* Symbols_sync(i : i+length(sync_sym)-1)));
    end

    figure; plot(corval); title(['Sync Correlation Peak (', mod_scheme, ')']);
    [~, start_sym_ind] = max(corval);

    % Phase Ambiguity & Demod
    preamble_rx = Symbols_sync(start_sym_ind : start_sym_ind + length(sync_sym) - 1);
    phase_ambiguity = angle(sum(preamble_rx .* conj(sync_sym)));
    Symbols_corrected = Symbols_sync * exp(-1i * phase_ambiguity);

    figure; scatter(real(Symbols_corrected), imag(Symbols_corrected), 'r.');
    title(['Constellation Diagram (', mod_scheme, ')']); grid on; axis square;

    num_symbols_expected = length(tx_ref_bits) / bps;
    rx_packet_sym = Symbols_corrected(start_sym_ind : start_sym_ind + num_symbols_expected - 1);
    bits_rec = qamdemod(rx_packet_sym, M, 'OutputType', 'bit', 'UnitAveragePower', true);
    
    % Align starting index for parsing below
    start_bit_ind = 1; 
end

%% 6. Extract Header and Image Data

if strcmp(mod_scheme, 'BPSK')
    idx = start_bit_ind + N;     % Dynamic BPSK offset
else
    idx = N + 1;                 % QAM bits_rec is already truncated exactly to packet start
end

% % --- Extract Rows (16 bits) ---
% rows_str = char(bits_rec(idx : idx + 15) + '0')';
% rx_rows  = bin2dec(rows_str);
% idx      = idx + 16;
% 
% % --- Extract Columns (16 bits) ---
% cols_str = char(bits_rec(idx : idx + 15) + '0')';
% rx_cols  = bin2dec(cols_str);
% idx      = idx + 16;
% ------------------------------------------------
% --- Extract Rows (bit-wise majority voting) ---
rows_rep = bits_rec(idx : idx + 48 - 1);

rows_majority = zeros(16,1);
for i = 1:16
    bits_triplet = rows_rep(3*i-2 : 3*i);   % take 3 repeated bits
    rows_majority(i) = sum(bits_triplet) >= 2;
end

rows_str = char(rows_majority + '0')';
rx_rows = bin2dec(rows_str);

idx = idx + 48;

% --- Extract Columns (bit-wise majority voting) ---
cols_rep = bits_rec(idx : idx + 48 - 1);

cols_majority = zeros(16,1);
for i = 1:16
    bits_triplet = cols_rep(3*i-2 : 3*i);
    cols_majority(i) = sum(bits_triplet) >= 2;
end

cols_str = char(cols_majority + '0')';
rx_cols = bin2dec(cols_str);

idx = idx + 48;

% --- Display Header ---
disp(['Received Image Header indicates size: ', ...
      num2str(rx_rows), 'x', num2str(rx_cols)]);

% --- Extract Image Data ---
img_data_length = rx_rows * rx_cols * 8;
rx_img_bits     = bits_rec(idx : idx + img_data_length - 1);%% 7. Bit Error Rate (BER) Calculation
if strcmp(mod_scheme, 'BPSK')
    total_packet_length = N + 48 + 48 + img_data_length; % Use N here
    received_packet = bits_rec(start_bit_ind : start_bit_ind + total_packet_length - 1);
else
    received_packet = bits_rec(1:length(tx_ref_bits));
end

csvwrite('received_packet.csv', received_packet); 

Nof_err_bits = sum(abs(tx_ref_bits - received_packet)); 
BER = Nof_err_bits / length(tx_ref_bits);

disp(['Total Bits Transmitted: ', num2str(length(tx_ref_bits))]);
disp(['Number of Error Bits: ', num2str(Nof_err_bits)]);
disp(['Bit Error Rate (BER): ', num2str(BER)]);
%% 8. Reconstruct and Display Image
img_bin_rec = reshape(rx_img_bits, 8, []).'; 
img_dec = bin2dec(char(img_bin_rec + '0'));

img_rec_matrix = reshape(img_dec, rx_rows, rx_cols);
final_image = uint8(img_rec_matrix);

figure;
imshow(final_image);
title(['Received Image | ', mod_scheme, ' | BER: ', num2str(BER)]);