clearvars -except 
close all

mes_num=10000;

for ii=1:10
    
 ii   
 res=zeros(mes_num,7);
 
 Depreciation_percentage=0.1;
p=1-exp(-Depreciation_percentage);
serviv=0.005+0.02*ii;
cell_num=1;


flur_start=1000;%round(sqrt(10)^ii)+3 ;
add_noise=false;
noise_amp=1;
intencity_for_singel_flur=753;

var_parameter=serviv;


    for jj=1:mes_num
        
 
         
        
frame_num=round(-log(serviv)/Depreciation_percentage);

%%%%%%%%%%%%%%%%%%%%%   build_flur_matrix
[flur,frame_num]=build_flur_matrix(frame_num,cell_num,flur_start,p,add_noise,noise_amp,intencity_for_singel_flur );
%%%%%%%%%%%%%%%%%%%%%%

% damp  = colcolate_damp( flur,cell_num,frame_num);
 
 damp=1-p;
 b_check(ii,jj)=damp+p-1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
meas_delta=ones(frame_num-1,cell_num);
expectet_delta=ones(frame_num-1,cell_num);
noise=ones(frame_num-1,cell_num);

for j=1:cell_num
    for i=1:frame_num-1
   
      meas_delta(i,j)=flur(i,j)-flur(i+1,j);
      expectet_delta(i,j)=flur(i,j)-flur(i,j)*damp;
      noise(i,j)= expectet_delta(i,j)- meas_delta(i,j);
     
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


noise2=noise;
noise=noise.^2;
 

if add_noise==true
  noise_factor=((noise(1:end,:)-flur(1:end-1,:)*(damp)-flur(2:end,:))./(expectet_delta(1:end,:)*(damp)));
else
  noise_factor=(noise(1:end,:)./(expectet_delta(1:end,:)*(damp)));  
end
 noise_anlaze=mean(noise_factor');
      noise_factor_mean=mean(noise_factor(1:end,:));
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  artical_res  = artical_new( flur,damp,add_noise );
  artical_res_without_integral  = artical_res_no_integral( flur,damp,add_noise );
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





res(jj,1)=mean(noise_factor_mean)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(jj,2)=var(noise_factor_mean)/intencity_for_singel_flur^2;
res(jj,3)=mean(artical_res)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(jj,4)=var(artical_res)/intencity_for_singel_flur^2;
res(jj,5)=var_parameter;
res(jj,6)=mean(artical_res_without_integral)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(jj,7)=var(artical_res_without_integral)/intencity_for_singel_flur^2;
    end

mean_v(ii)=mean(res(:,1));
mean_v2(ii)=mean(res(:,3));  
mean_v2_artical_res_no_integral(ii)=mean(res(:,6));
var_v(ii)=var(res(:,1));
var_v2(ii)=var(res(:,3)); 
var_v2_artical_res_no_integral(ii)=var(res(:,6)); 






v_m(ii)=var_parameter;    
    
end



%T = table(mean_rate_for_singel_flur,noise_factor_mean,noise_factor_abs,artical_res,result_fit_all_data,gause_fit,'VariableNames',{'rate_for_singel_flur_main' 'noise_factor_mean' 'noise_factor_abs' 'artical_res' 'result_fit_all_data' 'gause_fit'})

summery(1,1)=mean(noise_factor_mean);
summery(1,2)=mean(artical_res);


 summery(2,1)=var(noise_factor_mean);
 summery(2,2)=var(artical_res);



%table(summery','RowNames',{'rate_for_singel_flur_main' 'noise_factor_mean' 'noise_factor_abs' 'artical_res' 'result_fit_all_data' 'gause_fit'})
%table(summery','RowNames',{'noise_factor_mean'  'artical_res' })
%res

figure(2)
errorbar(100*v_m,mean_v,var_v)
hold on
errorbar(100*v_m,mean_v2,var_v2)
 errorbar(100*v_m,mean_v2_artical_res_no_integral,var_v2_artical_res_no_integral)     
legend('PFC','FFC','FFC with no integral')
xlabel('Percentage of survivors [%]')
ylabel('\nu^{2}/\nu_{est}')
mean(res)
