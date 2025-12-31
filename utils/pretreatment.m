function[num_frame,Ns,len_winframe,win] = pretreatment(xin,x,fs,M)
%% Framing & windowing
%setup for framing and windows
len_win = 0.016;       % length of window 16 ms
shift_percent = 0.25;       % win_shift 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% framing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len_winframe = fix(len_win * fs);
% window = ones(len_winframe,1);
window = hamming(len_winframe);
y = cell(M,1);
for i = 1:M  % generate input signals and noises for all nodes
    [y{i,1}, num_frame] = KFrame(x(:,i), len_winframe, window, shift_percent);
end

N = len_winframe;
K = len_winframe;
Ns = len_winframe*shift_percent;

win  = hamming(K);
win = win/sqrt(sum(win(1:Ns:K).^2));

% for j = 1:M    
%     X_win(:,:,j) = fft(enframe(x(:,j),win,Ns)');   
% end
frame_inf{1} = xin; frame_inf{2} = len_winframe; frame_inf{3} = shift_percent; frame_inf{4} = num_frame;
end