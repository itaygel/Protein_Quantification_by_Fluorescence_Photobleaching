function [ d ] = artical( flur,damp,add_noise )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
[frame_num,cell_num]=size(flur);
for i=1:cell_num
    for j=1:frame_num
        flur_idial(j,i)=flur(1,i)*damp^(j-1);
       

end
end
 flur_noise=flur_idial-flur;
 flur_noise=flur_noise.^2;
 
 for i=1:cell_num
    for j=1:frame_num
        flur_noise(j,i)=flur_noise(j,i)/(flur(1,i)*damp^(j-1)*(1-damp^(j-1)));
    end
 end
 
 d=mean(flur_noise(2:end,:));
     
 
 