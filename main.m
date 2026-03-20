clc
clear

addpath('input')
addpath('converters')
addpath('framing')
addpath('transmission')
addpath('utils')

% Select file
[file,path] = uigetfile({'*.mpeg;*.jpeg;*.txt'},'Select file');
filepath = fullfile(path,file);

% Load file and convert to binary
[binary_data,type] = load_file(filepath);

% Framing
payload = 320;
frames = frame_builder(binary_data,payload,type);

disp("Total Frames:")
disp(length(frames))

% Transmission stage
tx_signal = transmit_frames(frames);
disp("Transmission done")
