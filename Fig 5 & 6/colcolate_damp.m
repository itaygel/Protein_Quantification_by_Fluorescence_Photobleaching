function [ damp ] = colcolate_damp( flur,cell_num,frame_num)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
for j=1:cell_num
    for i=2:frame_num
   
      delta_calibrate(i,j)=flur(i,j)/flur(i-1,j);
      
     
    end
end
delta_calibrate=delta_calibrate(2:frame_num,:);
%damp=mean2(-log(delta_calibrate)/log(exp(1)));
damp=(mean(mean(delta_calibrate)));

end

