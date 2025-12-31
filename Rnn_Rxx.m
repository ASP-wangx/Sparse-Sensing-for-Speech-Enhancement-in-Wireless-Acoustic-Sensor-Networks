function [R,frame_idx] = Rnn_Rxx(input, frame_inf)

xin=frame_inf{1};len_frame= frame_inf{2}; overlap = frame_inf{3} ; Len = frame_inf{4} ; 
M = size(input,2);
x = input;
N = len_frame;
K = len_frame;
Ns = len_frame*overlap;
%% VAD
for i = 1:M
    x_vad(:,i) = awgn(xin(:,i),20,'measured');
    thh(i) = Vad_thr(x_vad(:,i),len_frame, overlap); % threshold
end

%%  Rvv and Rxx computation
flag = 0;
frame_idx = zeros(Len,1);

Rvv = zeros(M,M,K/2+1);
Ryy = zeros(M,M,K/2+1);
Rxx = zeros(M,M,K/2+1);

% Rvv 
for i = 1:Len
    
    for j = 1:M
        x_tmp = x_vad((i-1)*Ns+(1:N),j);
        P_tmp = 1/N * sum(x_tmp.^2);
        if P_tmp >= thh
            break;
        end
    end    
    
    if j == M  % all nodes capture noise frame  
        flag = flag + 1;
        frame_idx(i) = 1;
        
        for j = 1:M
            X(:,j) = fft(x((i-1)*Ns+(1:N),j),K);  % FFT
        end
                  
        for k = 1:K/2+1  
             corrvar = X(k,:).' * conj(X(k,:));
            Rvv(1:M,1:M,k) = Rvv(1:M,1:M,k) + corrvar/ trace(corrvar);            
        end 
    end      
end

for i = 1:Len   
    for j = 1:M
        X(:,j) = fft(x((i-1)*Ns+(1:N),j),K);  % FFT
    end
    for k = 1:K/2+1  
         corrvar = X(k,:).' * conj(X(k,:));
        Ryy(1:M,1:M,k) = Ryy(1:M,1:M,k) + corrvar/ trace(corrvar);            
    end 
end

Rvv = Rvv/flag;
Ryy = Ryy/Len;
Rxx = Ryy - Rvv;
R{1} = Ryy; R{2} = Rvv; R{3} = Rxx;
end