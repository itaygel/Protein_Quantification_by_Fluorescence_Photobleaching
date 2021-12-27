function [ flur,frame_num ] = build_flur_matrix(frame_num,cell_num,flur_start,p,noise,noise_amp,intencity_for_singel_flur,noise_bg_amp )
%UNTITLED3 Summary of this function goes here
%  בונה מערך של ערכי פלורוצניה כפי שנקראו בגלאי לפי מספר הפראמים ולפי מספר
%  החיידקים
%   Detailed explanation goes here
flur=zeros(frame_num,cell_num);
flur(1,:)=(flur(1,:)+1)*flur_start;%+(randn(1,cell_num))*0.05*flur_start;
%flur(1,:)=(flur(1,:)+1)*flur_start+(randn(1,cell_num))*0.05*flur_start;


for j=1:cell_num
    for i=1:frame_num-1
      bliching=0;  
      for u=1:flur(i,j)  
if rand<=p
    bliching=bliching+1;
end
      end
     flur(i+1,j) =flur(i,j) -  bliching; 
    % bliching_sum(i,j)=bliching;
    %%%%%%%%%%%%%%%%%%%%%%%
    if flur(i+1,j)==0
        frame_num=i;
        break
    end
    
    
    
    
    end
end

flur=flur*intencity_for_singel_flur;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%add noise 
if noise==true
    bg_amp=[];
    for i=1:frame_num
  bg_amp=[flur(1,:);bg_amp];  
    end
    
X = randn(size(flur));
rr=X.*sqrt(flur);
flur=flur+rr*noise_amp;

X_bg = randn(size(flur)).*sqrt(bg_amp*noise_bg_amp);
flur=flur+X_bg;
end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fit

end

