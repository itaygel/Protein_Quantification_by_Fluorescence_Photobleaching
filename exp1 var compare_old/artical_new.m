function [ v ] = artical_new( flur,damp,add_noise )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
[frame_num,cell_num]=size(flur);

p_end=1-damp^(frame_num-1);

for i=1:cell_num
    for j=1:frame_num
        flur_idial(j,i)=flur(1,i)*damp^(j-1);
       

end
end
 flur_noise=flur_idial-flur;
 %
 %flur_noise=flur_noise.^2-flur(1,1);
 if add_noise==true
 flur_noise=flur_noise.^2-2*flur_idial;
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
vv(j-1)=flur_noise(j,i)*((1-damp^(j-1+0)-(1-damp^(j-2+0)))/2+(1-damp^(j+0)-(1-damp^(j-1+0)))/2);
  
        
    end
   v(i)=sum(vv)/(flur(1,i)*((3*p_end^2-2*p_end^3)/6)); 

end




     
 
 