clearvars -except a
close all
mes_num=1;

for iii=1:1
for ii=1:5
    res=zeros(mes_num,7);
    
 ii   ;
 
tdtu=0.05;
p=1-exp(-tdtu);
serviv=0.1;
cell_num=50;
flur_start=1000; 
add_noise=false;
noise_amp=1;
intencity_for_singel_flur=753*tdtu/0.1;

var_parameter=flur_start;
 
    for jj=1:mes_num

frame_num=round(-log(serviv)/tdtu);




%%%%%%%%%%%%%%%%%%%%%
[flur,frame_num]=build_flur_matrix(frame_num,cell_num,flur_start,p,add_noise,noise_amp,intencity_for_singel_flur );

frame_numa(jj,ii)=frame_num;
%%%%%%%%%%%%%%%%%%%%%%

 
 
 damp=(1-p)*(0.998+0.002*ii);

 b_check(ii,jj)=damp+p-1;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 var_parameter1(ii)=0.999+0.001*ii; 
 
 
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
   noise_factor2=(noise2(1:end,:)./(expectet_delta(1:end,:)*(damp)));  

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fit line to all data
if add_noise==true
  noise_factor=((noise(1:end,:)-flur(1:end-1,:)*(damp)-flur(2:end,:))./(expectet_delta(1:end,:)*(damp)));
else
  noise_factor=(noise(1:end,:)./(expectet_delta(1:end,:)*(damp)));  
end
 noise_anlaze=mean(noise_factor');
      noise_factor_mean=mean(noise_factor(1:end,:));
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
artical_res  = artical_new( flur,damp,add_noise );
artical_res_no_integral  = artical( flur,damp );
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot histogram





res(jj,1)=mean(noise_factor_mean)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(jj,2)=var(noise_factor_mean)/intencity_for_singel_flur^2;
res(jj,3)=mean(artical_res)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(jj,4)=var(artical_res)/intencity_for_singel_flur^2;
res(jj,5)=var_parameter;
res(jj,6)=mean(artical_res_no_integral)/intencity_for_singel_flur;%/intencity_for_singel_flur;
res(jj,7)=var(artical_res_no_integral)/intencity_for_singel_flur^2;
    end

  v_var_m(ii,iii)=mean(res(:,1));
v_var_m2(ii,iii)=mean(res(:,3));  
 v_var_m2_artical_res_no_integral(ii)=mean(res(:,6));
mean_var_m(ii)=var(res(:,1));
mean_var_m2(ii)=var(res(:,3)); 
mean_var_m2_artical_res_no_integral(ii)=var(res(:,6)); 

v_m(ii)=var_parameter;    
  p_error(ii,iii)=((flur_start*(1-p-damp)^2)/(damp-damp^2))*((1-p)^(frame_num-1)-1)/(-p*(frame_num-1))+(p-p^2)/(damp-damp^2);
  p_integral=1/(((1-serviv)^2)/2-((1-serviv)^3)/3);
 pp_error(ii,iii)=1+(2*flur_start*(0.002*(ii-1))^2);
end

end

summery(1,1)=mean(noise_factor_mean);
summery(1,2)=mean(artical_res);
summery(2,1)=var(noise_factor_mean);
summery(2,2)=var(artical_res);



figure(4)

hold on
plot(var_parameter1,v_var_m,'o')
plot(var_parameter1,v_var_m2,'*')
plot(var_parameter1,p_error)
xlabel('p_{error}')
ylabel('\nu/\nu_{real}')
legend('PFC','FFC','Error from equation')