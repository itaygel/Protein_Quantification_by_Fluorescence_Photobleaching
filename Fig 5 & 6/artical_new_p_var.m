function [ v ] = artical_new_p_var( flur,dam,add_noise,noise_amp,noise_bg_amp )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
[frame_num,cell_num]=size(flur);
%damp=mean(flur')/mean(flur(1,:));
%damp=(flur')/mean(flur(1,:));
%damp(1)=[];
%p_end=1-damp^(frame_num-1);

%damp=mean(flur')/mean(flur(1,:));
for k=1:frame_num
    d(k)=dam.^(k-1);
end
damp=d;
p_end=1-damp(end);



for i=1:cell_num
    for j=1:frame_num
        flur_idial(j,i)=flur(1,i)*damp(j);
       

    end
end
 flur_noise=flur_idial-flur;
 %
 %flur_noise=flur_noise.^2-flur(1,1);
 if add_noise==true
 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  bg noise fix
%      ooo=ones(1,cell_num);
%     oo=[];
%     for i1=1:frame_num
%   oo=[oo;ooo*damp^i1];  
%     end
     bg_amp=[];
    for i=1:frame_num
  bg_amp=[bg_amp;flur(1,:)];  
    end
    bg_amp=bg_amp*noise_bg_amp;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
     %fix  2/05/2018
 flur_noise=flur_noise.^2-2*flur_idial*noise_amp-bg_amp-flur_idial*noise_bg_amp;

 
 
 else
   flur_noise=flur_noise.^2; 
 end
 %flur_noise=flur_noise.^2;
 
 %%%%%%%%%%%%%%%%%%%%%% trapz
% for i=1:cell_num
%     
%     for j=1:frame_num
%  y(j)=flur_noise(j,i)/((damp)^(j-1));
%         x(j)=1-(damp)^(j-1);
%     end
%     v(i)=trapz(y,x)/flur(1,i);
% end
%  
% s=mean(flur_noise');
% plot(x,s)

 %%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 

for i=1:cell_num
    
    for j=2:frame_num
 
%vv(j-1)=flur_noise(j,i)*(1-damp^(j-1+0)-(1-damp^(j-2+0)));
%vv(j-1)=flur_noise(j,i)*((1-damp^(j-1+0)-(1-damp^(j-2+0)))/2+(1-damp^(j+0)-(1-damp^(j-1+0)))/2);
if j==frame_num
    vv(j-1)=flur_noise(j,i)*((1-damp(j)-(1-damp(j-1)))/2+(1-damp(j)-(1-damp(j)))/2);
    

else
vv(j-1)=flur_noise(j,i)*((1-damp(j)-(1-damp(j-1)))/2+(1-damp(j+1)-(1-damp(j)))/2);
end
        
    end
   v(i)=sum(vv)/(flur(1,i)*((3*p_end^2-2*p_end^3)/6)); 

end




     
 
 