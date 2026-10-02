% ==========================================
% --- TRANSMITTER: Image over RF (Hybrid) ---
% ==========================================
clear; clc;

%% --- USER CONFIGURATION ---
% Create a pop-up dialog for Modulation Scheme
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

%% 1. Read and Process Image
img = imread('/MATLAB Drive/new codes/lowkbinput image.jpeg'); 
if size(img, 3) == 3
    img = rgb2gray(img); 
end

img = imresize(img,[64 64]);

rows = size(img, 1);
cols = size(img, 2);

img_bin = dec2bin(img(:), 8);
img_bin_t = img_bin.';
img_bits = img_bin_t(:) - '0';

%% 2. Create Header (Sync Bits + Rows + Cols)
h_pn = comm.PNSequence('Polynomial', [6 5 0], 'InitialConditions', [0 0 0 0 0 1], 'VariableSizeOutput', true, 'MaximumOutputSize', [N, 1]);
sync_bits = h_pn(N);

rows_bin = dec2bin(rows, 16).';
rows_bits = rows_bin(:) - '0';

cols_bin = dec2bin(cols, 16).';
cols_bits = cols_bin(:) - '0';

%% 3. Build the Packet
% tx_packet = [sync_bits; rows_bits; cols_bits; img_bits];
% -------- my code--------------
rows_bits_rep = reshape(repmat(rows_bits.', 3, 1), [], 1);
cols_bits_rep = reshape(repmat(cols_bits.', 3, 1), [], 1);

tx_packet = [sync_bits; rows_bits_rep; cols_bits_rep; img_bits];

%% 4. Modulation & Pulse Shaping
switch filter_type
    case 'Rect',  FIR_coeff = ones(samplesPerSymbol, 1) / samplesPerSymbol;
    case 'RRC25', FIR_coeff = rcosdesign(0.25, 8, samplesPerSymbol, 'sqrt').';
    case 'RRC50', FIR_coeff = rcosdesign(0.50, 8, samplesPerSymbol, 'sqrt').';
    case 'RRC75', FIR_coeff = rcosdesign(0.75, 8, samplesPerSymbol, 'sqrt').';
end

if strcmp(mod_scheme, 'BPSK')
    % --- ORIGINAL BPSK LOGIC ---
    bits_I = 2 * (tx_packet - 0.5); 
    sig_ipt_I = upsample(bits_I, samplesPerSymbol);
    
    tx_signal = conv(FIR_coeff, sig_ipt_I);
    tx_signal = tx_signal / max(abs(tx_signal));
    
    csvwrite('Transmitted_ref.csv', tx_packet); 
    csvwrite('Transmitted_data.csv', tx_signal); % 1D array for LabVIEW
    
else
    % --- QAM LOGIC ---
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
    csvwrite('Transmitted_data.csv', [real(tx_signal), imag(tx_signal)]); % 2D array (I & Q)
end

disp(['--- TX Complete ---']);
disp(['Modulation: ', mod_scheme, ' | Filter: ', filter_type]);
disp(['Image Size: ', num2str(rows), 'x', num2str(cols)]);