function z = dr_mvdr_ss(input, frame_inf, sour, F_msc, cons_inf, slc, A, Rvv, frame_idx)

%% input
 xin=frame_inf{1};len_frame= frame_inf{2}; overlap = frame_inf{3} ; Len = frame_inf{4} ; 
[att,dly]=steer(sour);
att = att(slc);
dly = dly(slc);
glo_W  = cons_inf{1}; T0 =  cons_inf{2}; 

A1 = A(slc,slc);
W = EEW_opt_symm(A1);

M = size(input,2);

x = input;
xs = x(:,slc);

N = len_frame;
K = len_frame;
Ns = len_frame*overlap;

Rvv = Rvv(slc,slc,:);

R = zeros(size(xs,2),size(xs,2),K/2+1); IR0 = zeros(size(xs,2),size(xs,2),K/2+1);
for k = 1:K/2+1  
    R(:,:,k) = Rvv(:,:,k);
    IR0(:,:,k) = inv( R(:,:,k));
end
%% window function
win  = hamming(K);
win = win/sqrt(sum(win(1:Ns:K).^2));
%% MVDR beamforming
for j = 1:size(xs,2)    
    X_win(:,:,j) = fft(enframe(xs(:,j),win,Ns)');   
end

alpha = 0.990; 
XR = zeros(K,size(xs,2)); IR = zeros(size(xs,2),size(xs,2),K/2+1);
for i = 1:Len
    
    % FFT of received signal
    for j = 1:size(xs,2)
        XR(:,j) = fft(xs((i-1)*Ns+(1:N),j),K);  % FFT
    end
    
    if(frame_idx(i)==1)
        for k = 1:K/2+1 
            % atf
            dvec = att.*exp(-sqrt(-1)*2*pi*(k-1)/K*F_msc*dly); 
            % compute q0H by distributed consensus strategy
            temp_q = XR(k,:)' .* IR0(:,:,k);    
            for t = 1:500
                temp_q = W*temp_q;
            end
            q(:,k) = size(xs,2)*temp_q(1,:)';
            
            YR = inv(IR0(:,:,k))*q(:,k);  
              
           IR(:,:,k) = ( 1/alpha*IR0(:,:,k) - (  1/alpha*(1-alpha)*q(:,k)*q(:,k)' )/( alpha + (1-alpha)*q(:,k)'* YR )  );
            
            Xout(k,i) = ( dvec'*IR(:,:,k)*YR  )/( dvec'*IR(:,:,k)*dvec );
        end
        IR0 = IR;
    else
        for k = 1:K/2+1 
            % atf
            dvec = att.*exp(-sqrt(-1)*2*pi*(k-1)/K*F_msc*dly); 
            % compute q0H by distributed consensus strategy
            temp_q = XR(k,:)' .* IR0(:,:,k);        
            for t = 1:T0
                temp_q = W*temp_q;
            end
            q(:,k) = size(xs,2)*temp_q(1,:)';

           IR(:,:,k) = IR0(:,:,k);
            
           Xout(k,i) = ( dvec'*q(:,k) )/( dvec'*IR(:,:,k)*dvec );     
        end
    end 
    %% 2nd half frequencies
    for k = 2:K/2
        Xout(k+K/2,i) =  conj( Xout(K/2+2-k,i));
    end
end
%% output

xtmp = real(ifft(Xout))';
z = overlapadd(xtmp,win,Ns);
end