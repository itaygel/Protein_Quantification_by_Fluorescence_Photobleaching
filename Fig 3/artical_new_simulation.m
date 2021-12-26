function [ v ] = artical_new_simulation( damp,frame_num )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
cell_num=1;
p_end=damp^(frame_num-1);
flur(1,1)=1;
for i=1:cell_num
    for j=1:frame_num
        flur_idial(j,i)=flur(1,1)*damp^(j-1);
       

    end
end
 flur_noise=flur_idial;

   flur_noise=flur_noise; 
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
   v(i)=sum(vv)/(flur(1,1)*((3*p_end^2-2*p_end^3)/6)); 
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   p_start=0.905;
  

   1/((1.5*p_end^2-1.5*p_start^2)/(p_end^3-p_start^3)-1)
   
   

end




     
 
 