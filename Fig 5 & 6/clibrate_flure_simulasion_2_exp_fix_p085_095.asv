clearvars -except 
close all
mes_num=10;
fit_for_all=0

fix_noise=false;%true; ;
add_noise=false;%true;%true;
noise_amp=0;%*ii-0.1;
%noise_bg_amp=0.5;

noise_bg_amp=0.0;
num_of_test=1;


% fix_noise=true; 
% add_noise=true;%true;
% noise_amp=1;%*ii-0.1;
% %noise_bg_amp=0.5;
% 
% noise_bg_amp=0.1;
% num_of_test=1;


pf1=0.75;
for ii=1:10
    
    res=zeros(mes_num,7);
    ii
    
 
tdtu=0.05;
p=1-exp(-tdtu);
serviv=0.01 ;
cell_num=30;
%flur_start=1000; %round(sqrt(10)^7) ;

intencity_for_singel_flur=753;%25*ii;%*ii;%*(tdtu/0.1)*ii/10;

var_parameter=intencity_for_singel_flur;
 
    for jj=1:mes_num

%frame_num=round(-log(serviv)/tdtu);%10;%*qqq;




%%%%%%%%%%%%%%%%%%%%%
 %[flur,frame_num]=build_flur_matrix(frame_num,cell_num,flur_start,p,add_noise,noise_amp,intencity_for_singel_flur,noise_bg_amp );
% 
% frame_numa(jj,ii)=frame_num;
%%%%%%%%%%%%%%%%%%%%%%
nnn(ii)=round(500*10^(ii/5));
nnn(ii)=5000;

[ flur,frame_num ] = build_flur_matrix_p_var(15+ii*3,cell_num,nnn(ii),1,false,0,intencity_for_singel_flur,0 );
%p_var(ii)=1-0.99+ii*0.02;
p_var(ii)=1-pf1-0.1+ii*0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ yt1,frame_num ] = build_flur_matrix(16,cell_num,round(0.5*nnn(ii)),p_var(ii),false,0,750,0 );
% %%[ yt1,frame_num ] = build_flur_matrix(15,cell_num,round(0.5*nnn(ii)),1-0.7,false,0,750,0 );
% [ yt2,frame_num ] = build_flur_matrix(16,cell_num,nnn(ii)-round(0.5*nnn(ii)),1-pf1,false,0,750,0 );
% %[ yt3,frame_num ] = build_flur_matrix(15,cell_num,nnn(ii)-round(0.5*nnn(ii)),0.1,false,0,750,0 );
% 
% %flur=yt2;
% flur=yt1+yt2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %damp  = colcolate_damp( flur,cell_num,frame_num);

   
   
 %damp_pfc=mean(flur(2:end)').\mean(flur(1:end-1)');
 
 
 if fit_for_all==1
        [fitresult, gof] =two_exp_fit(mean(flur'));
 
 
         %ey=feval(fitresult, 1:frame_num);
        ey=feval(fitresult, 1:frame_num)-fitresult.f;
        %ey=ey*y1/10;
        
        flur=flur-fitresult.f;

   
   damp=ey/ey(1);
   damp_pfc=ey(2:end)./ey(1:end-1);
   end
   
   
   
   
 %damp=(1-p)*(0.999+0.001*ii);   for p error
 %damp=(1-p);
 %b_check(ii,jj)=damp+p-1;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %var_parameter1(ii)=0.999+0.001*ii; 
 var_parameter1(ii)=ii;
 
 meas_delta=ones(frame_num-1,cell_num);
expectet_delta=ones(frame_num-1,cell_num);
noise=ones(frame_num-1,cell_num);

for j=1:cell_num
    
    
     if fit_for_all==0
        [fitresult, gof] =two_exp_fit(flur(:,j));
 
 
         %ey=feval(fitresult, 1:frame_num);
        ey=feval(fitresult, 1:frame_num)-fitresult.f;
        %ey=ey*y1/10;
        
        flur=flur-fitresult.f;

   
   damp=ey/ey(1);
   damp_pfc=ey(2:end)./ey(1:end-1);
   end
    
    
    
    
    
    for i=1:frame_num-1
   
      meas_delta(i,j)=flur(i,j)-flur(i+1,j);
      %%%%expectet_delta=(flur(i,j)-exp_fit(j,3))*(1-exp(-exp_fit(j,2)))
      %expectet_delta(i,j)=(flur(i,j)-exp_fit(j,3))-(flur(i,j)-exp_fit(j,3))*(exp(-exp_fit(j,2)));
      %expectet_delta(i,j)=(flur(i,j)-exp_fit_mean_c(qqq,iiii))-(flur(i,j)-exp_fit_mean_c(qqq,iiii))*(exp(-exp_fit_mean_b(qqq,iiii)));
      %expectet_delta(i,j)=(flur(i,j)-exp_fit(j,3))-(flur(i,j)-exp_fit(j,3))*(exp(-exp_fit_mean_b(qqq,iiii)));
      expectet_delta(i,j)=flur(i,j)-flur(i,j)*damp_pfc(i);%(exp(-damp));
      noise(i,j)= expectet_delta(i,j)- meas_delta(i,j);
     
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


noise2=noise;
 noise=noise.^2;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fit line to all data
% bin_num=40;
% % [result_fit_all_data_t_p,result_fit_all_data_t_p_b, a1,b1]= liner_fit_to_noise( expectet_delta,noise,bin_num,damp );
% result_fit_all_data(ii,jj)=1;%result_fit_all_data_t_p;
% result_fit_all_data_b(ii,jj)=1;%result_fit_all_data_t_p_b;
% result_fit_a(ii,jj)=1;%a1;
% result_fit_b(ii,jj)=1;%b1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%noise=noise-noise*result_fit_all_data_t_p_b/mean2(noise);
%noise=mean2(noise)- (result_fit_all_data_t_p_b);
 
 
 %noise_factor_abs(ii,jj)=1.57*mean2((abs(noise2(1:end-1,:))./sqrt(expectet_delta(2:end,:)))).^2;
%  noise_factor_abs(ii,jj)=1.57*mean2((abs(noise2(1:end,:))./sqrt(expectet_delta(1:end,:)*(damp)))).^2;    
%  
%  fix_noise=(noise(1:end,:)-flur(2:end,:));
%  fix_noise(fix_noise<0)=0;
%  noise_factor=(fix_noise./(expectet_delta(1:end,:)*(damp)));
 


if fix_noise==true
    
    bg_amp=[];
    for i=1:frame_num-1
  bg_amp=[bg_amp;flur(1,:)];  
    end
    bg_amp=bg_amp*noise_bg_amp;
    
 % noise_factor=((noise(1:end,:)-flur(1:end-1,:)*(damp)-flur(2:end,:)-bg_amp)./(expectet_delta(1:end,:)*(damp)));
% noise_factor=((noise(1:end,:)-2*flur(1:end-1,:)*(damp)-bg_amp)./(expectet_delta(1:end,:)*(damp)));

%noise_factor=((noise(1:end,:)-bg_amp-bg_amp*damp)./(expectet_delta(1:end,:)*(damp)))-2/(1-damp);
noise_factor=((noise(1:end,:)-bg_amp-bg_amp*damp-2*flur(1:end-1,:)*(damp)*noise_amp)./(expectet_delta(1:end,:)*(damp)));

%noise_factor=((noise(1:end,:))./(expectet_delta(1:end,:)*(damp)))-(2/(1-damp))*noise_amp-noise_bg_amp*(1+damp)/(1-damp);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%תיקון מסדר שני
% fix2=2*flur(2:end,:)./((((noise_factor*0.5-flur(2:end,:))/damp)+1).*(damp*(1-damp)));
% noise_factor=noise_factor+fix2;

%noise_factor=noise_factor-(noise_bg_amp/((1-damp))).*(1+1./oo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
  noise_factor=(noise(1:end,:)./(expectet_delta(1:end,:).*(damp_pfc)));  
end
 noise_anlaze=mean(noise_factor');
   %   noise_factor_mean(ii,jj)=mean2(noise_factor(1:end,:));
      noise_factor_mean=mean(noise_factor(1:end,:));
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
artical_res  = artical_new_p_var( flur,damp,fix_noise,noise_amp,noise_bg_amp );
 artical_res_no_integral  = artical_2_exp( flur,damp,fix_noise,noise_amp,noise_bg_amp );
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot histogram
% %noise2_factor=(noise2(1:end-1,:)./expectet_delta(2:end,:));
% %noise2_factor=(noise2(1:end,:)./expectet_delta(1:end,:));
% noise2_factor=(noise2(1:end,:)./sqrt(expectet_delta(1:end,:)));
% noise2_factor=noise2_factor*sqrt(mean(expectet_delta(1,:)));
% 
% %noise2_factor=noise2_factor*mean(expectet_delta(1,:));
% % 
% % [N,edges] = histcounts(noise2_factor);
% % edges=edges(1:end-1);
% % edges=edges+edges(2)/2-edges(1)/2;
% % %% Fit: 'untitled fit 1'.
% % [xData, yData] = prepareCurveData( edges, N );
% % 
% % % Set up fittype and options.
% % ft = fittype( 'gauss1' );
% % opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% % opts.Display = 'Off';
% % opts.Lower = [-Inf -Inf 0];
% % opts.Robust = 'Bisquare';
% % opts.StartPoint = [317 -10000 39240.6496472761];
% % 
% % 
% % % excludedPoints = excludedata( xData, yData, 'Indices', [46 47] );
% % % opts.Exclude = excludedPoints;
% % 
% % 
% % 
% % % Fit model to data.
% % [fitresult, gof] = fit( xData, yData, ft, opts );
% % % 
% % % 
% %  fit_a(ii)=fitresult.a1;
% %  fit_b(ii)=fitresult.b1;
% %  fit_c(ii)=fitresult.c1/2;
%  %bbb(ii)=fit_c(ii)/((flur(1,1))*p*(1-p))
%  gause_fit(ii)=1;%4.0323*fit_c(ii)^2/((flur(1,1))*damp*(1-damp)*2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




res(jj,1)=mean(noise_factor_mean)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(jj,2)=var(noise_factor_mean)/intencity_for_singel_flur^2;
res(jj,3)=mean(artical_res)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(jj,4)=var(artical_res)/intencity_for_singel_flur^2;
res(jj,5)=var_parameter;
res(jj,6)=mean(artical_res_no_integral)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(jj,7)=var(artical_res_no_integral)/intencity_for_singel_flur^2;
    end

  v_var_m(ii)=mean(res(:,1));
v_var_m2(ii)=mean(res(:,3));  
 v_var_m2_artical_res_no_integral(ii)=mean(res(:,6));
mean_var_m(ii)=var(res(:,1));
mean_var_m2(ii)=var(res(:,3)); 
mean_var_m2_artical_res_no_integral(ii)=var(res(:,6)); 

v_m(ii)=var_parameter;    
 %p_error(ii)=((flur_start*(1-p-damp)^2)/(damp-damp^2))*((1-p)^(frame_num)-1)/(-p*frame_num)+(p-p^2)/(damp-damp^2);
%  p_error(ii)=((flur_start*(1-p-damp)^2)/(damp-damp^2))*((1-p)^(frame_num-1)-1)/(-p*(frame_num-1))+(p-p^2)/(damp-damp^2);

   
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
% % 
% % R = normrnd(0,5,10000,1);
% % [N,edges] = histcounts(R);
% %  edges=edges(1:end-1);
% % edges=edges+edges(2)/2-edges(1)/2;
% figure(2)
% errorbar(v_m,v_var_m,mean_var_m)
% hold on
% errorbar(v_m,v_var_m2,mean_var_m2)
%       
% 
% 
% mean(res)
% figure(3)
% semilogx(v_m,mean_var_m2,'*b')
% hold on
% %plot(v_m,0.8-3./(5*v_m),'*g')
% semilogx(v_m,0.8-3./(5*v_m),'*g')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot for p error
% figure(4)
% plot(var_parameter1,p_error)
% hold on
% plot(var_parameter1,v_var_m,'o')
% plot(var_parameter1,v_var_m2,'*')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% save(['ex2_fix_' 'noise_bg_amp_' num2str(noise_bg_amp*100) '_noise_amp' num2str(noise_amp*100)])
% 
% 
% clearvars -except noise_amp noise_bg_amp mes_num num_of_test
% close all
% %mes_num=10000;
% 
% 
% fix_noise=true; %false;
% add_noise=true;%false%true;


% 
% 
% for ii=1:num_of_test
%   
%     res=zeros(mes_num,7);
%     ii
%     
%  
% tdtu=0.05;
% p=1-exp(-tdtu);
% serviv=0.1 ;
% cell_num=50;
% flur_start=1000; %round(sqrt(10)^7) ;
% 
% intencity_for_singel_flur=25*ii;%*ii;%*(tdtu/0.1)*ii/10;
% 
% var_parameter=intencity_for_singel_flur;
%  
%     for jj=1:mes_num
% 
% frame_num=round(-log(serviv)/tdtu);%10;%*qqq;
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%
% [flur,frame_num]=build_flur_matrix(frame_num,cell_num,flur_start,p,add_noise,noise_amp,intencity_for_singel_flur,noise_bg_amp );
% 
% frame_numa(jj,ii)=frame_num;
% %%%%%%%%%%%%%%%%%%%%%%
% 
%  damp  = colcolate_damp( flur,cell_num,frame_num);
%  
%  %damp=(1-p)*(0.999+0.001*ii);   for p error
%  %damp=(1-p);
%  b_check(ii,jj)=damp+p-1;
% 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %var_parameter1(ii)=0.999+0.001*ii; 
%  
%  
%  meas_delta=ones(frame_num-1,cell_num);
% expectet_delta=ones(frame_num-1,cell_num);
% noise=ones(frame_num-1,cell_num);
% 
% for j=1:cell_num
%     for i=1:frame_num-1
%    
%       meas_delta(i,j)=flur(i,j)-flur(i+1,j);
%       %%%%expectet_delta=(flur(i,j)-exp_fit(j,3))*(1-exp(-exp_fit(j,2)))
%       %expectet_delta(i,j)=(flur(i,j)-exp_fit(j,3))-(flur(i,j)-exp_fit(j,3))*(exp(-exp_fit(j,2)));
%       %expectet_delta(i,j)=(flur(i,j)-exp_fit_mean_c(qqq,iiii))-(flur(i,j)-exp_fit_mean_c(qqq,iiii))*(exp(-exp_fit_mean_b(qqq,iiii)));
%       %expectet_delta(i,j)=(flur(i,j)-exp_fit(j,3))-(flur(i,j)-exp_fit(j,3))*(exp(-exp_fit_mean_b(qqq,iiii)));
%       expectet_delta(i,j)=flur(i,j)-flur(i,j)*damp;%(exp(-damp));
%       noise(i,j)= expectet_delta(i,j)- meas_delta(i,j);
%      
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% noise2=noise;
%  noise=noise.^2;
%  
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fit line to all data
% % bin_num=40;
% % % [result_fit_all_data_t_p,result_fit_all_data_t_p_b, a1,b1]= liner_fit_to_noise( expectet_delta,noise,bin_num,damp );
% % result_fit_all_data(ii,jj)=1;%result_fit_all_data_t_p;
% % result_fit_all_data_b(ii,jj)=1;%result_fit_all_data_t_p_b;
% % result_fit_a(ii,jj)=1;%a1;
% % result_fit_b(ii,jj)=1;%b1;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %noise=noise-noise*result_fit_all_data_t_p_b/mean2(noise);
% %noise=mean2(noise)- (result_fit_all_data_t_p_b);
%  
%  
%  %noise_factor_abs(ii,jj)=1.57*mean2((abs(noise2(1:end-1,:))./sqrt(expectet_delta(2:end,:)))).^2;
% %  noise_factor_abs(ii,jj)=1.57*mean2((abs(noise2(1:end,:))./sqrt(expectet_delta(1:end,:)*(damp)))).^2;    
% %  
% %  fix_noise=(noise(1:end,:)-flur(2:end,:));
% %  fix_noise(fix_noise<0)=0;
% %  noise_factor=(fix_noise./(expectet_delta(1:end,:)*(damp)));
%  
% 
% 
% if fix_noise==true
%     
%     bg_amp=[];
%     for i=1:frame_num-1
%   bg_amp=[bg_amp;flur(1,:)];  
%     end
%     bg_amp=bg_amp*noise_bg_amp;
%     
%  % noise_factor=((noise(1:end,:)-flur(1:end-1,:)*(damp)-flur(2:end,:)-bg_amp)./(expectet_delta(1:end,:)*(damp)));
% % noise_factor=((noise(1:end,:)-2*flur(1:end-1,:)*(damp)-bg_amp)./(expectet_delta(1:end,:)*(damp)));
% 
% noise_factor=((noise(1:end,:)-bg_amp-bg_amp*damp-2*flur(1:end-1,:)*(damp)*noise_amp)./(expectet_delta(1:end,:)*(damp)));
% %noise_factor=((noise(1:end,:)-bg_amp-bg_amp*damp)./(expectet_delta(1:end,:)*(damp)))-2/(2-damp);
% 
% 
% % noise_factor=(noise(1:end,:)./(expectet_delta(1:end,:)*(damp)))-2/(1-damp);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%תיקון מסדר שני
% % fix2=2*flur(2:end,:)./((((noise_factor*0.5-flur(2:end,:))/damp)+1).*(damp*(1-damp)));
% % noise_factor=noise_factor+fix2;
% 
% %noise_factor=noise_factor-(noise_bg_amp/((1-damp))).*(1+1./oo);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% else
%   noise_factor=(noise(1:end,:)./(expectet_delta(1:end,:)*(damp)));  
% end
%  noise_anlaze=mean(noise_factor');
%    %   noise_factor_mean(ii,jj)=mean2(noise_factor(1:end,:));
%       noise_factor_mean=mean(noise_factor(1:end,:));
%       
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% artical_res  = artical_new( flur,damp,fix_noise,noise_amp,noise_bg_amp );
%  artical_res_no_integral  = artical( flur,damp,fix_noise,noise_amp,noise_bg_amp );
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% res(jj,1)=mean(noise_factor_mean)/intencity_for_singel_flur;%/intencity_for_singel_flur;
% res(jj,2)=var(noise_factor_mean)/intencity_for_singel_flur^2;
% res(jj,3)=mean(artical_res)/intencity_for_singel_flur;%/intencity_for_singel_flur;
% res(jj,4)=var(artical_res)/intencity_for_singel_flur^2;
% res(jj,5)=var_parameter;
% res(jj,6)=mean(artical_res_no_integral)/intencity_for_singel_flur;%/intencity_for_singel_flur;
% res(jj,7)=var(artical_res_no_integral)/intencity_for_singel_flur^2;
%     end
% 
%   v_var_m(ii)=mean(res(:,1));
% v_var_m2(ii)=mean(res(:,3));  
%  v_var_m2_artical_res_no_integral(ii)=mean(res(:,6));
% mean_var_m(ii)=var(res(:,1));
% mean_var_m2(ii)=var(res(:,3)); 
% mean_var_m2_artical_res_no_integral(ii)=var(res(:,6)); 
% 
% v_m(ii)=var_parameter;    
%  %p_error(ii)=((flur_start*(1-p-damp)^2)/(damp-damp^2))*((1-p)^(frame_num)-1)/(-p*frame_num)+(p-p^2)/(damp-damp^2);
%   p_error(ii)=((flur_start*(1-p-damp)^2)/(damp-damp^2))*((1-p)^(frame_num-1)-1)/(-p*(frame_num-1))+(p-p^2)/(damp-damp^2);
% 
%    
% %   corect_noise=2/(1-damp);
% % v_fix=mean(noise_factor_mean)+corect_noise*(1-mean(noise_factor_mean)/(mean(noise_factor_mean)+corect_noise))
% %  
% 
% end
% % % figure(1)
% % % errorbar(res(:,5),res(:,1),res(:,2))
% % % hold on
% % % errorbar(res(:,5),res(:,3),res(:,4))
% % 
% % % result
% % %gause_fit=gause_fit';
% % %T = table(mean_rate_for_singel_flur,noise_factor_mean,noise_factor_abs,artical_res,result_fit_all_data,gause_fit,'VariableNames',{'rate_for_singel_flur_main' 'noise_factor_mean' 'noise_factor_abs' 'artical_res' 'result_fit_all_data' 'gause_fit'})
% % %summery(1,1)=mean(mean_rate_for_singel_flur);
% % summery(1,1)=mean(noise_factor_mean);
% % %summery(1,3)=mean(noise_factor_abs);
% % summery(1,2)=mean(artical_res);
% % % summery(1,5)=mean(result_fit_all_data);
% % % summery(1,6)=mean(gause_fit);
% % 
% % 
% % % summery(2,1)=var(mean_rate_for_singel_flur);
% %  summery(2,1)=var(noise_factor_mean);
% % % summery(2,3)=var(noise_factor_abs);
% %  summery(2,2)=var(artical_res);
% % % summery(2,5)=var(result_fit_all_data);
% % % summery(2,6)=var(gause_fit);
% % 
% % 
% % %table(summery','RowNames',{'rate_for_singel_flur_main' 'noise_factor_mean' 'noise_factor_abs' 'artical_res' 'result_fit_all_data' 'gause_fit'})
% % %table(summery','RowNames',{'noise_factor_mean'  'artical_res' })
% % %res
% % % 
% % % R = normrnd(0,5,10000,1);
% % % [N,edges] = histcounts(R);
% % %  edges=edges(1:end-1);
% % % edges=edges+edges(2)/2-edges(1)/2;
% % figure(2)
% % errorbar(v_m,v_var_m,mean_var_m)
% % hold on
% % errorbar(v_m,v_var_m2,mean_var_m2)
% %       
% % 
% % 
% % mean(res)
% % figure(3)
% % semilogx(v_m,mean_var_m2,'*b')
% % hold on
% % %plot(v_m,0.8-3./(5*v_m),'*g')
% % semilogx(v_m,0.8-3./(5*v_m),'*g')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot for p error
% figure(4)
% plot(var_parameter1,p_error)
% hold on
% plot(var_parameter1,v_var_m,'o')
% plot(var_parameter1,v_var_m2,'*')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  for ex 2
% plot(v_m,v_var_m)
% hold on
% plot(v_m,v_var_m2)
% 
% 
% 
% % Create xlabel
% xlabel({'\nu'});
% 
% % Create ylabel
% ylabel({'\nu_{error}/\nu_{real}'});
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% op_damp=1/damp;
% bg_error=(noise_bg_amp*(1+damp)/(1-damp))*(op_damp^frame_num-1)/((op_damp-1)*frame_num)
% int_limit=1/((1-serviv)^2/2-(1-serviv)^3/3);
% bg_error1=((1-serviv)+(1-serviv)^2/2)*int_limit*noise_bg_amp;
% %bg_error1=(1.5-(1-serviv)-(1-serviv)^2/2)*int_limit*noise_bg_amp;
% 
% plot(v_m,(v_m+bg_error)./v_m,'k')
% plot(v_m,(v_m+bg_error1)./v_m,'k')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% save(['ex2_no_fix_' 'noise_bg_amp_' num2str(noise_bg_amp*100) '_noise_amp' num2str(noise_amp*100)])

p_var2=(1-p_var)/pf1;

plot(p_var2,v_var_m,'o')
hold on
plot(nnn,v_var_m2_artical_res_no_integral,'*')
errorbar(nnn,v_var_m,sqrt(mean_var_m))
errorbar(nnn,v_var_m2,sqrt(mean_var_m2_artical_res_no_integral))
errorbar(nnn,v_var_m2_artical_res_no_integral,sqrt(mean_var_m2_artical_res_no_integral))


xlabel('p_2/p_1')
ylabel('\nu_{est}/\nu')

figure()
hold on
plot(nnn,v_var_m,'o')
plot(nnn,v_var_m2)
plot(nnn,v_var_m2_artical_res_no_integral)
save ([ 'Two_exp_' num2str(fit_for_all), '.mat'])