function [x,xin,isnr,av_pow3,fs] = add_noise(RIR,RIR1,RIR2,SNR,SIR1,SIR2,len)

M = size(RIR,1);
[signal,fs]=audioread('./speech/speech_female.wav');  %clean speech 
[noise1,fs]=audioread('./speech/noise1.wav');  %interference 1 

N = len*fs;

noise2 = randn(N,1);     % interference 2
signal=signal(1:N); noise1=noise1(1:N);



av_pow = mean( sum(signal(:,1).^2,1)/size(signal(:,1),1) );  
av_pow1 = mean( sum(noise1(:,1).^2,1)/size(noise1(:,1),1) );  
av_pow2 = mean( sum(noise2(:,1).^2,1)/size(noise2(:,1),1) );  

sigma_n1 = sqrt( av_pow/av_pow1/(10^(SIR1/10)) );
sigma_n2 = sqrt( av_pow/av_pow2/(10^(SIR2/10)) );	



for i = 1:M  % generate input signals and noises for all nodes
    xin(:,i) = filter(RIR{i,1}(:,1),1,signal);
    ns1(:,i) = filter(RIR1{i,1}(:,1),1,sigma_n1*noise1);
    ns2(:,i) = filter(RIR2{i,1}(:,1),1,sigma_n2*noise2);
end
mic_noise = randn(N,1);     %ambient noise 
av_pow = mean( sum(xin(:,1).^2,1)/size(xin(:,1),1) );  
av_pow3 = mean( sum(mic_noise(:,1).^2,1)/size(mic_noise(:,1),1) );  

sigma = sqrt( av_pow/av_pow3/(10^(SNR/10)) );	

for i = 1:M    
    noise(:,i) =sigma* randn(N,1)+ ns1(:,i) + ns2(:,i);
    s_xin(i) = (xin(:,i)'*xin(:,i))/length(N);
    s_noise(i) = (noise(:,i)'*noise(:,i))/length(N);
    iSNR(i) = 10*log10( s_xin(i)/s_noise(i) );
    x(:,i) = xin(:,i) +  noise(:,i);     
end
isnr = mean(iSNR);



end

