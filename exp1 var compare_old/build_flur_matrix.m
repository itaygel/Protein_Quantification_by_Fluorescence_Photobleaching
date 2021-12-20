function [ flur,frame_num ] = build_flur_matrix(frame_num,cell_num,flur_start,p,noise,noise_amp,intencity_for_singel_flur )
%UNTITLED3 Summary of this function goes here

% Constructs an array of fluorescence values as read in the detector
% The size of the array is the number of frames * the number of bacteria


flur=zeros(frame_num,cell_num);
flur(1,:)=(flur(1,:)+1)*flur_start;%+(randn(1,cell_num))*0.05*flur_start;




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
        frame_num=i+1;
        break
    end
    
    
    
    
    end
end

flur=flur(1:frame_num)*intencity_for_singel_flur;
%flur=flur*intencity_for_singel_flur;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%add noise 
if noise==true
X = randn(size(flur));
rr=X.*sqrt(flur);
flur=flur+rr*noise_amp;
end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fit

end

