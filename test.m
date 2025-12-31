clc
clear all

%% setup
load('./inputs/mic_position.mat'); % load mic positions
M = size(mic,1);   % mic number
addpath("./spares_ss/");
addpath("./utils/");


%% Optimal symmetric
thr = 2.05;
A = Graph_gene(thr,mic);
W = EEW_opt_symm(A);
A0 = ceil(abs(W));
A0 = A0 - diag(diag(A0));   % Adjacency matrix

%% parameters 
p = 15;
b0 = 0.4;
eta0 = 0.1;
delta  = 1;
belta = 3.5;
K0 = 128; 
Q = 1; 
Num = 3;
SNR = 20; % dB  (ambient noise)
SIR1 = 5; % dB (interference 1)
SIR2 = 5; % dB (interference 2)
len = 5;        % signal length (s)

%% source positions
sour = [2.5,6];
itf1 = [2,3];
itf2 = [7.5,2];
sources{1} = sour; sources{2} = itf1; sources{3} = itf2;

load("./inputs/rt200.mat");
RIR = RIR_cell;
load("./inputs/rt200noise1.mat");
RIR1 = RIR_cell;
load("./inputs/rt200noise2.mat");
RIR2 = RIR_cell;

[x,xin,isnr,av_pow3,fs] = add_noise(RIR,RIR1,RIR2,SNR,SIR1,SIR2,len);
N = len*fs;


%% Rnn_Rxx computing
T0 = ceil(log(1/1e-4) /log(1/0.7595));
cons_inf{1} = W; cons_inf{2} = T0;
len_win = 0.016;       % length of window 10 ms
shift_percent = 0.25;       % win_shift 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% framing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len_winframe = fix(len_win * fs);
window = rectwin(len_winframe);

x0=x(1:N,:); xin0 = xin(1:N,:);
[num_frame,Ns,K,win] = pretreatment(xin0,x0,fs,M);

frame_inf{1} = xin0; frame_inf{2} = K; frame_inf{3} = Ns/K; frame_inf{4} = num_frame; 
[R,frame_idx] = Rnn_Rxx(x0, frame_inf);

x2 = x(1:len*fs,:);
xin2 = xin(1:len*fs,:); 
[yy, num_frame2] = KFrame(x2, len_winframe, window, shift_percent);%重叠分帧
frame_inf2{1} = xin2; frame_inf2{2} = K; frame_inf2{3} = Ns/K; frame_inf2{4} = num_frame2;

[R2,frame_idx2] = Rnn_Rxx(x2, frame_inf2);

%% proposed D-FI-SS alg
[d_i,d_w,d_choose,d_t,lambda2] = SS_drib_12inf_R(A, R, sources, b0, K0, Q,belta);

d_index = find(d_choose>=eta0);

%% distributed mvdr
d_out = dr_mvdr_ss(x2, frame_inf2, sources{1}, fs, cons_inf, d_index, A, R2{2}, frame_idx2);






