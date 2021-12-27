function [ new_yt ] = zip_yt( yt,n )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

[e,l]=size(yt);
new_frame_num=floor(l/n);

for i=1:new_frame_num
   
    new_yt(1:e,i)=sum(yt(1:e,(i-1)*n+1:(i)*n),2);
    
    
    
    
    
    
end



end

