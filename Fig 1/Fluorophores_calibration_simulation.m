clearvars -except 
close all
mes_num=1;


fix_noise=false; %false;
add_noise=false;%false%true;
noise_amp=1;%*ii-0.1;
%noise_bg_amp=0.5;

noise_bg_amp=0.25;
num_of_test=15;
 res=zeros(mes_num,7);
for ii=1:num_of_test
    
   
    ii
    
 
tdtu=0.1;
p=1-exp(-tdtu);
serviv=0.1;%0.003+(ii-1)*0.005 ;
cell_num=10000;
flur_start=1000; %round(sqrt(10)^7) ;
flur_srart_var=ii*0.01;
intencity_for_singel_flur=753;%*ii;%*(tdtu/0.1)*ii/10;

var_parameter=intencity_for_singel_flur;
 
    for jj=1:mes_num

frame_num=round(-log(serviv)/tdtu);%10;%*qqq;




%%%%%%%%%%%%%%%%%%%%%
[flur,frame_num]=build_flur_matrix(frame_num,cell_num,flur_start,p,add_noise,noise_amp,intencity_for_singel_flur,noise_bg_amp,flur_srart_var );

frame_numa(jj,ii)=frame_num;
%%%%%%%%%%%%%%%%%%%%%%

 damp  = colcolate_damp( flur,cell_num,frame_num);
 

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
 

if fix_noise==true
    
    bg_amp=[];
    for i=1:frame_num-1
  bg_amp=[bg_amp;flur(1,:)];  
    end
    bg_amp=bg_amp*noise_bg_amp;
    

noise_factor=((noise(1:end,:)-bg_amp-bg_amp*damp-2*flur(1:end-1,:)*(damp)*noise_amp)./(expectet_delta(1:end,:)*(damp)));





%%%%%%%%%%%%%%
% fix2=2*flur(2:end,:)./((((noise_factor*0.5-flur(2:end,:))/damp)+1).*(damp*(1-damp)));
% noise_factor=noise_factor+fix2;

%noise_factor=noise_factor-(noise_bg_amp/((1-damp))).*(1+1./oo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
  noise_factor=(noise(1:end,:)./(expectet_delta(1:end,:)*(damp)));  
end
 noise_anlaze=mean(noise_factor');
      noise_factor_mean=mean(noise_factor(1:end,:));
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
artical_res  = artical_new( flur,damp,fix_noise,noise_amp,noise_bg_amp );
artical_res_no_integral  = artical( flur,damp,fix_noise,noise_amp,noise_bg_amp );
artical_res_no_p  = artical_new_new( flur,damp,fix_noise,noise_amp,noise_bg_amp );
artical_res_no_p_fix  = artical_new_new_fix( flur,damp,fix_noise,noise_amp,noise_bg_amp );

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





res(ii,1)=mean(noise_factor_mean)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(ii,2)=var(noise_factor_mean)/intencity_for_singel_flur^2;
res(ii,3)=mean(artical_res)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(ii,4)=var(artical_res)/intencity_for_singel_flur^2;
res(ii,5)=var_parameter;
res(ii,6)=mean(artical_res_no_integral)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(ii,7)=var(artical_res_no_integral)/intencity_for_singel_flur^2;
res(ii,8)=mean(artical_res_no_p)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(ii,9)=var(artical_res_no_p)/intencity_for_singel_flur^2;   
res(ii,10)=mean(artical_res_no_p_fix)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(ii,11)=var(artical_res_no_p_fix)/intencity_for_singel_flur^2;   
    end

%   v_var_m(ii)=mean(res(:,1));
% v_var_m2(ii)=mean(res(:,3));  
%  v_var_m2_artical_res_no_integral(ii)=mean(res(:,6));
%  
% mean_var_m(ii)=var(res(:,1));
% mean_var_m2(ii)=var(res(:,3)); 
% mean_var_m2_artical_res_no_integral(ii)=var(res(:,6)); 
% mean_var_m2_artical_res_no_p(ii)=var(res(:,8)); 
% v_m(ii)=var_parameter;    
 %p_error(ii)=((flur_start*(1-p-damp)^2)/(damp-damp^2))*((1-p)^(frame_num)-1)/(-p*frame_num)+(p-p^2)/(damp-damp^2);
  p_error(ii)=((flur_start*(1-p-damp)^2)/(damp-damp^2))*((1-p)^(frame_num-1)-1)/(-p*(frame_num-1))+(p-p^2)/(damp-damp^2);

   
%   corect_noise=2/(1-damp);
% v_fix=mean(noise_factor_mean)+corect_noise*(1-mean(noise_factor_mean)/(mean(noise_factor_mean)+corect_noise))
%  

end
% % figure(1)
% % errorbar(res(:,5),res(:,1),res(:,2))
% % hold on
% % errorbar(res(:,5),res(:,3),res(:,4))
% 
% % result
% %gause_fit=gause_fit';
% %T = table(mean_rate_for_singel_flur,noise_factor_mean,noise_factor_abs,artical_res,result_fit_all_data,gause_fit,'VariableNames',{'rate_for_singel_flur_main' 'noise_factor_mean' 'noise_factor_abs' 'artical_res' 'result_fit_all_data' 'gause_fit'})
% %summery(1,1)=mean(mean_rate_for_singel_flur);
% summery(1,1)=mean(noise_factor_mean);
% %summery(1,3)=mean(noise_factor_abs);
% summery(1,2)=mean(artical_res);
% % summery(1,5)=mean(result_fit_all_data);
% % summery(1,6)=mean(gause_fit);
% 
% 
% % summery(2,1)=var(mean_rate_for_singel_flur);
%  summery(2,1)=var(noise_factor_mean);
% % summery(2,3)=var(noise_factor_abs);
%  summery(2,2)=var(artical_res);
% % summery(2,5)=var(result_fit_all_data);
% % summery(2,6)=var(gause_fit);
% 
% 
% %table(summery','RowNames',{'rate_for_singel_flur_main' 'noise_factor_mean' 'noise_factor_abs' 'artical_res' 'result_fit_all_data' 'gause_fit'})
% %table(summery','RowNames',{'noise_factor_mean'  'artical_res' })
% %res

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot for p error
% figure(4)
% plot(var_parameter1,p_error)
% hold on
% plot(var_parameter1,v_var_m,'o')
% plot(var_parameter1,v_var_m2,'*')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(['ex2_fix_' 'noise_bg_amp_' num2str(noise_bg_amp*100) '_noise_amp' num2str(noise_amp*100)])

res
close all
figure(1)

set(0,'DefaultFigureWindowStyle' , 'normal')
set(gcf,'PaperUnits','cent','PaperPosition',[0 0 7.1 6.5])

e=0.01:0.01:ii*0.01;
plot(e,res(:,8),'*')
hold on
p_start=(1+damp)/2;
p_end=damp^frame_num;%serviv;

%plot(e,(e./(1-e)).^2*flur_start*intencity_for_singel_flur*(-p_start-log(1-p_start)+p_end+log(1-p_end))/intencity_for_singel_flur+1,':')
%plot(e,(e).^2*flur_start*intencity_for_singel_flur*(-p_start-log(1-p_start)+p_end+log(1-p_end))/intencity_for_singel_flur+1,':')
%plot(e,(e).^2*flur_start*intencity_for_singel_flur*1.44/intencity_for_singel_flur+1,':')
%((p_start^3)/3-(p_end^3)/3)/((p_start^2)/2-(p_start^3)/3-(p_end^2)/2+(p_end^3)/3)
res(:,8)'./(e.^2*flur_start*intencity_for_singel_flur*(-p_start-log(1-p_start)+p_end+log(1-p_end))/(intencity_for_singel_flur)+1)
plot(e,res(:,10),'*')
plot(e,(e).^2*flur_start*intencity_for_singel_flur*((p_start^3)/3-(p_end^3)/3)/((p_start^2)/2-(p_start^3)/3-(p_end^2)/2+(p_end^3)/3)/intencity_for_singel_flur+1,':')

p_start=0.9;%%%%%%%%%%%%%  test

%plot(e,(e).^2*flur_start*intencity_for_singel_flur*(p_end-p_start+log(1-p_end)-log(1-p_start))/intencity_for_singel_flur+1,':')
%%%%%%%%%%%%%%%%%%
ylabel({'\nu_{error}/\nu_{real}'});
xlabel({'\alpha'})
legend('Flur var','P var')
%print(['my_figure1'], '-dtiff', '-r900')