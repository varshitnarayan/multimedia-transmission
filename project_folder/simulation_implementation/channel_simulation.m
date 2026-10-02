% ==========================================
% --- THE VIRTUAL CABLE (Lab Simulator) ---
% Run this AFTER the Transmitter, BEFORE the Receiver
% ==========================================
% clear; clc;
disp('Simulating transmission through the air...');

% 1. Read what the hybrid transmitter made
tx_signal = csvread('Transmitted_data.csv');

% 2. Check if it's BPSK (1 col) or QAM (2 cols)
if size(tx_signal, 2) == 1
    I_channel = tx_signal(:,1); 
    Q_channel = zeros(length(I_channel), 1); % BPSK has no Q channel natively
else
    I_channel = tx_signal(:,1);
    Q_channel = tx_signal(:,2);
end

% --- to add channel noise ---
% % This will force your receiver's PLL to actually do some work!
% I_channel = I_channel + 0.05 * randn(size(I_channel));
% Q_channel = Q_channel + 0.05 * randn(size(Q_channel));
% % ------------------------------------

% 3. Combine them into two columns (Col 1 = I, Col 2 = Q)
simulated_rx_data = [I_channel, Q_channel];

% 4. Save it as the file the Receiver is looking for
csvwrite('Acquired_data.csv', simulated_rx_data);

disp('Lab simulated successfully! You can now run your Receiver code.');