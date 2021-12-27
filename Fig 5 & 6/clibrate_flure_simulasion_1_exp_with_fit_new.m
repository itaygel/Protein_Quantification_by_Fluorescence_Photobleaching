clearvars -except flur
close all

mes_num=500;
fit_for_all=0

fix_noise=false;%true; ;
add_noise=false;%true;%true;
noise_amp=0;%*ii-0.1;
noise_bg_amp=0.0;
num_of_test=1;


for ii=1:20 %35
    
    res=zeros(mes_num,11);
    ii
    
 p=0.15;
serviv=0.1 ;
cell_num=1;
intencity_for_singel_flur=753;%25*ii;%*ii;%*(tdtu/0.1)*ii/10;
var_parameter=intencity_for_singel_flur;
 
    for jj=1:mes_num




%%%%%%%%%%%%%%%%%%%%%
 %[flur,frame_num]=build_flur_matrix(frame_num,cell_num,flur_start,p,add_noise,noise_amp,intencity_for_singel_flur,noise_bg_amp );
%
%frame_num=32
% frame_numa(jj,ii)=frame_num;
%%%%%%%%%%%%%%%%%%%%%%
%nnn(ii)=round(120*10^(ii/5));
nnn(ii)=5000;

[ flur,frame_num ] = build_flur_matrix(3+(ii-1)*2,cell_num,nnn(ii),p,false,0,intencity_for_singel_flur,0 );
%[ flur,frame_num ] = build_flur_matrix(3,cell_num,nnn(ii),p,false,0,intencity_for_singel_flur,0 );
%[ flur,frame_num ] = build_flur_matrix(10,cell_num,nnn(ii),p+(ii-1)/20,false,0,intencity_for_singel_flur,0 );

%p_var(ii)=1-0.99+ii*0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %damp  = colcolate_damp( flur,cell_num,frame_num);

   
   
 %damp_pfc=mean(flur(2:end)').\mean(flur(1:end-1)');
 
 
 if fit_for_all==1
        [fitresult, gof] =one_exp_fit(mean(flur'));
        ey=feval(fitresult, 1:frame_num)-fitresult.f;        
        flur=flur-fitresult.f;

   
   damp=ey/ey(1);
   damp_pfc=ey(2:end)./ey(1:end-1);
   end
   
   
   
 
meas_delta=ones(frame_num-1,cell_num);
expectet_delta=ones(frame_num-1,cell_num);
noise=ones(frame_num-1,cell_num);

for j=1:cell_num
    
     if fit_for_all==0
        [fitresult, gof] =one_exp_fit(flur(:,j));
 p_from_fit(jj,ii)=fitresult.b;
 
         %ey=feval(fitresult, 1:frame_num);
        ey=feval(fitresult, 1:frame_num)-fitresult.f;
        %ey=ey*y1/10;
        
        flur=flur-fitresult.f;

   
   damp=ey/ey(1);
   damp_pfc=ey(2:end)./ey(1:end-1);
   end
    
    
    
    
    
    for i=1:frame_num-1
%    
%       meas_delta(i,j)=flur(i,j)-flur(i+1,j);
       expectet_delta(i,j)=flur(i,j)-flur(i,j)*damp_pfc(i);
      noise(i,j)= flur(i+1,j)-flur(i,j)*damp_pfc(i);
     
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% expected from fit only
if fit_for_all==0
noise1=(flur(1:end,j)-ey).^2;
%noise=noise(2:end,:);
pp=(mean(damp_pfc));
%noise_factor2=noise(1:end,:)./(flur(1:end-1,:)*(pp*(1-pp)));  

for k=1:frame_num-1
   noise_factor2(k,j)=noise1(k+1)./(flur(1,j)*(pp^k*(1-pp^k)));  

    
end
%noise_factor2=noise(1:end,:)./(flur(:,1:end-1)*(pp*(1-pp)));  
for k=1:frame_num-1
   noise_factor_pfc(k,j)=noise1(k+1)./(flur(k,j)*(pp*(1-pp)));  

    
end




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    






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
    
 % noise_factor=((noise(1:end,:)-flur(1:end-1,:)*(damp)-flur(2:end,:)-bg_amp)./(expectet_delta(1:end,:)*(damp)));
% noise_factor=((noise(1:end,:)-2*flur(1:end-1,:)*(damp)-bg_amp)./(expectet_delta(1:end,:)*(damp)));

%noise_factor=((noise(1:end,:)-bg_amp-bg_amp*damp)./(expectet_delta(1:end,:)*(damp)))-2/(1-damp);
noise_factor=((noise(1:end,:)-bg_amp-bg_amp*damp-2*flur(1:end-1,:)*(damp)*noise_amp)./(expectet_delta(1:end,:)*(damp)));

%noise_factor=((noise(1:end,:))./(expectet_delta(1:end,:)*(damp)))-(2/(1-damp))*noise_amp-noise_bg_amp*(1+damp)/(1-damp);





else
  noise_factor=(noise(1:end,:)./(expectet_delta(1:end,:).*(damp_pfc)));  
end
 noise_anlaze=mean(noise_factor');
   %   noise_factor_mean(ii,jj)=mean2(noise_factor(1:end,:));
      noise_factor_mean=mean(noise_factor(1:end,:));
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
artical_res  = artical_new_p_var( flur,mean(damp_pfc),fix_noise,noise_amp,noise_bg_amp );
%artical_res_no_integral  = artical_2_exp( flur,mean(damp_pfc),fix_noise,noise_amp,noise_bg_amp );
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% expected from fit only
if fit_for_all==1
noise=(flur(1:end,:)-ey).^2;
noise=noise(2:end,:);
pp=(mean(damp_pfc));
%noise_factor2=noise(1:end,:)./(flur(1:end-1,:)*(pp*(1-pp)));  

for k=1:frame_num-1
   noise_factor2(k,1:cell_num)=noise(k,:)./(flur(1,:)*(pp^k*(1-pp^k)));  

    
end
%noise_factor2=noise(1:end,:)./(flur(:,1:end-1)*(pp*(1-pp)));  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 



res(jj,1)=mean(noise_factor_mean)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(jj,2)=var(noise_factor_mean)/intencity_for_singel_flur^2;
res(jj,3)=mean(artical_res)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(jj,4)=var(artical_res)/intencity_for_singel_flur^2;
res(jj,5)=var_parameter;
%res(jj,6)=mean(artical_res_no_integral)/intencity_for_singel_flur;%/intencity_for_singel_flur;
%res(jj,7)=var(artical_res_no_integral)/intencity_for_singel_flur^2;
  res(jj,8)=mean((noise_factor2(:)))/intencity_for_singel_flur;
res(jj,9)=var(noise_factor2)/intencity_for_singel_flur^2;
res(jj,10)=mean(noise_factor_pfc)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(jj,11)=var(noise_factor_pfc)/intencity_for_singel_flur^2;

  
    
    
    end

  v_var_m(ii)=mean(res(:,1));
v_var_m2(ii)=mean(res(:,3));  
 v_var_m2_artical_res_no_integral(ii)=mean(res(:,6));
mean_var_m(ii)=var(res(:,1));
mean_var_m2(ii)=var(res(:,3)); 
mean_var_m2_artical_res_no_integral(ii)=var(res(:,6)); 
v_from_fit(ii)=mean(res(:,8));
v_from_fit_pfc(ii)=mean(res(:,10));


v_m(ii)=var_parameter;    
 %p_error(ii)=((flur_start*(1-p-damp)^2)/(damp-damp^2))*((1-p)^(frame_num)-1)/(-p*frame_num)+(p-p^2)/(damp-damp^2);
%  p_error(ii)=((flur_start*(1-p-damp)^2)/(damp-damp^2))*((1-p)^(frame_num-1)-1)/(-p*(frame_num-1))+(p-p^2)/(damp-damp^2);

   
%   corect_noise=2/(1-damp);
% v_fix=mean(noise_factor_mean)+corect_noise*(1-mean(noise_factor_mean)/(mean(noise_factor_mean)+corect_noise))
%  

end

%plot(p_var2,v_var_m,'o')
hold on
plot(nnn,v_var_m2_artical_res_no_integral,'*')
errorbar(nnn,v_var_m,sqrt(mean_var_m))
errorbar(nnn,v_var_m2,sqrt(mean_var_m2_artical_res_no_integral))
errorbar(nnn,v_var_m2_artical_res_no_integral,sqrt(mean_var_m2_artical_res_no_integral))


xlabel('p_2/p_1')
ylabel('\nu_{est}/\nu')

figure()
hold on
plot(v_var_m,'o')
plot(v_var_m2)
plot(v_var_m2_artical_res_no_integral)
plot(v_from_fit)
plot(v_from_fit_pfc)

legend('PFC', 'FFC' ,'FFC without integral' ,'FFC from fit','pfc from fit')
xlabel('Frame #')
ylabel('\nu_{est}/\nu')

save ([ 'one_exp_' num2str(fit_for_all), '.mat'])

figure()
plot(sqrt(var(p_from_fit))/p_from_fit(1,1))

m=mean(p_from_fit);
plot(mean(abs(p_from_fit-m))./(1-m))