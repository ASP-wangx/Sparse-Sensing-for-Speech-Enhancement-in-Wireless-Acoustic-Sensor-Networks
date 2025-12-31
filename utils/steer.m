function [att,dly]=steer(sour)

s = sour;


load('mic_position.mat');
M = size(mic,1);

 for i = 1:M
     att(i,1) = 1/norm(s-mic(i,:),2);
     dly(i,1) = norm(s-mic(i,:),2)/343;     
 end
%  att = att/att(10);
%  dly = dly-dly(10);
 
   
end