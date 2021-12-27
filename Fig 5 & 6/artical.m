function [ d ] = artical( flur,damp,add_noise,noise_amp,noise_bg_amp )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
[frame_num,cell_num]=size(flur);

%damp=mean(flur')/mean(flur(1,:));

for i=1:cell_num
    for j=1:frame_num
        flur_idial(j,i)=flur(1,i)*damp^(j-1);
       

end
end
 flur_noise=flur_idial-flur;
 %flur_noise=flur_noise.^2;
 
 
 if add_noise==true
 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  bg noise fix
     ooo=ones(1,cell_num);
    oo=[];
    for i1=1:frame_num
  oo=[oo;ooo*damp^i1];  
    end
     bg_amp=[];
    for i=1:frame_num
  bg_amp=[bg_amp;flur(1,:)];  
    end
    bg_amp=bg_amp*noise_bg_amp;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
     
 %flur_noise=flur_noise.^2-2*flur_idial-ones(size(flur_idial))*mean(flur(1,:));
 flur_noise=flur_noise.^2-2*flur_idial*noise_amp*damp-bg_amp--bg_amp*damp;
 else
   flur_noise=flur_noise.^2; 
 end
 
 
 
 
 for i=1:cell_num
    for j=1:frame_num
        flur_noise(j,i)=flur_noise(j,i)/(flur(1,i)*damp^(j-1)*(1-damp^(j-1)));
    end
 end
 
 d=mean(flur_noise(2:end,:));
     
 
 
 