% This program performs a simulation at first without noise correction,
% saves the data and then repeats the process with noise correction

clearvars -except 
close all
mes_num=200;


fix_noise=false;%true; ;
add_noise=true;%true;%true;
noise_amp=1; % Shot noise amp
%noise_bg_amp=0.5;

noise_bg_amp=0.25;
num_of_test=10;

for ii=1:num_of_test
    
    res=zeros(mes_num,7);
    ii
    
 
tdtu=0.05;
p=1-exp(-tdtu);
serviv=0.1 ;
cell_num=50;
flur_start=1000; 

intencity_for_singel_flur=75+25*ii;

var_parameter=intencity_for_singel_flur;
 
    for jj=1:mes_num

frame_num=round(-log(serviv)/tdtu);%10;%*qqq;




%%%%%%%%%%%%%%%%%%%%%
[flur,frame_num]=build_flur_matrix(frame_num,cell_num,flur_start,p,add_noise,noise_amp,intencity_for_singel_flur,noise_bg_amp );

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
  
      expectet_delta(i,j)=flur(i,j)-flur(i,j)*damp;%(exp(-damp));
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



else
  noise_factor=(noise(1:end,:)./(expectet_delta(1:end,:)*(damp)));  
end
 noise_anlaze=mean(noise_factor');
      noise_factor_mean=mean(noise_factor(1:end,:));
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
artical_res  = artical_new( flur,damp,fix_noise,noise_amp,noise_bg_amp );
 artical_res_no_integral  = artical( flur,damp,fix_noise,noise_amp,noise_bg_amp );
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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
  p_error(ii)=((flur_start*(1-p-damp)^2)/(damp-damp^2))*((1-p)^(frame_num-1)-1)/(-p*(frame_num-1))+(p-p^2)/(damp-damp^2);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(['ex2_fix_' 'noise_bg_amp_' num2str(noise_bg_amp*100) '_noise_amp' num2str(noise_amp*100)])


subplot(2,1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  for ex 2
plot(v_m,v_var_m)
hold on
plot(v_m,v_var_m2)



% Create xlabel
xlabel({'\nu'});

% Create ylabel
ylabel({'\nu_{error}/\nu_{real}'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%op_damp=1/damp;
%bg_error=(noise_bg_amp*(1+damp)/(1-damp))*(op_damp^frame_num-1)/((op_damp-1)*frame_num);
%int_limit=1/((1-serviv)^2/2-(1-serviv)^3/3);
%bg_error1=((1-serviv)+(1-serviv)^2/2)*int_limit*noise_bg_amp;
%bg_error1=(1.5-(1-serviv)-(1-serviv)^2/2)*int_limit*noise_bg_amp;

%plot(v_m,(v_m+bg_error)./v_m,'k')
%plot(v_m,(v_m+bg_error1)./v_m,'k:')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




int_start=damp;
int_end=serviv;

bg_error1=-noise_bg_amp*log((int_end)/(int_start))+2*noise_bg_amp*log((1-int_end)/(1-int_start))+2*log((int_end-1)/(int_start-1));


plot(v_m,(v_m+bg_error1)./v_m,'r:','LineWidth',1)

op_damp=1/damp;
bg_error=2/(1-damp)+(noise_bg_amp*(1+damp)/(1-damp))*(op_damp^frame_num-1)/((op_damp-1)*frame_num);

plot(v_m,(v_m+bg_error)./v_m,'k:','LineWidth',1)

set(gca, 'YGrid', 'on', 'XGrid', 'on')





legend('PFC','FFC','PFC Expected noise','FFC Expected noise')







clearvars -except noise_amp noise_bg_amp mes_num num_of_test



fix_noise=true; 
add_noise=true;




for ii=1:num_of_test
  
    res=zeros(mes_num,7);
    ii
    
 
tdtu=0.05;
p=1-exp(-tdtu);
serviv=0.1 ;
cell_num=50;
flur_start=1000; %round(sqrt(10)^7) ;

intencity_for_singel_flur=75+25*ii;%*ii;%*(tdtu/0.1)*ii/10;

var_parameter=intencity_for_singel_flur;
 
    for jj=1:mes_num

frame_num=round(-log(serviv)/tdtu);%10;%*qqq;




%%%%%%%%%%%%%%%%%%%%%
[flur,frame_num]=build_flur_matrix(frame_num,cell_num,flur_start,p,add_noise,noise_amp,intencity_for_singel_flur,noise_bg_amp );

frame_numa(jj,ii)=frame_num;
%%%%%%%%%%%%%%%%%%%%%%

 damp  = colcolate_damp( flur,cell_num,frame_num);
 
 %damp=(1-p)*(0.999+0.001*ii);   for p error
 %damp=(1-p);
 b_check(ii,jj)=damp+p-1;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %var_parameter1(ii)=0.999+0.001*ii; 
 
 
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
%noise_factor=((noise(1:end,:)-bg_amp-bg_amp*damp)./(expectet_delta(1:end,:)*(damp)))-2/(2-damp);



else
  noise_factor=(noise(1:end,:)./(expectet_delta(1:end,:)*(damp)));  
end
 noise_anlaze=mean(noise_factor');
      noise_factor_mean=mean(noise_factor(1:end,:));
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
artical_res  = artical_new( flur,damp,fix_noise,noise_amp,noise_bg_amp );
 artical_res_no_integral  = artical( flur,damp,fix_noise,noise_amp,noise_bg_amp );
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





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
  p_error(ii)=((flur_start*(1-p-damp)^2)/(damp-damp^2))*((1-p)^(frame_num-1)-1)/(-p*(frame_num-1))+(p-p^2)/(damp-damp^2);

   


end

subplot(2,1,2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  for ex 2
plot(v_m,v_var_m)
hold on
plot(v_m,v_var_m2)



% Create xlabel
xlabel({'\nu'});

% Create ylabel
ylabel({'\nu_{error}/\nu_{real}'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

op_damp=1/damp;
bg_error=(noise_bg_amp*(1+damp)/(1-damp))*(op_damp^frame_num-1)/((op_damp-1)*frame_num)
int_limit=1/((1-serviv)^2/2-(1-serviv)^3/3);
bg_error1=((1-serviv)+(1-serviv)^2/2)*int_limit*noise_bg_amp;
%bg_error1=(1.5-(1-serviv)-(1-serviv)^2/2)*int_limit*noise_bg_amp;

%plot(v_m,(v_m+bg_error)./v_m,'k')
%plot(v_m,(v_m+bg_error1)./v_m,'k:')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%legend('PFC','FFC','PFC Expected noise','FFC Expected noise')
save(['ex2_no_fix_' 'noise_bg_amp_' num2str(noise_bg_amp*100) '_noise_amp' num2str(noise_amp*100)])