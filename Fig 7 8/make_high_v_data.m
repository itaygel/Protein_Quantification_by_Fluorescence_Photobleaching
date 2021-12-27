for i=1:6
    
   yt=peter_Swain_data(i).yt;
f_name=peter_Swain_data(i).name;
for j=1:5
ytt(1:length(yt),j)=sum(yt(:,1+8*(j-1):1+8*(j))')';
end  
    
  
peter_Swain_data(i+19).yt=ytt;
peter_Swain_data(i+19).name=[f_name,'_8'];


end