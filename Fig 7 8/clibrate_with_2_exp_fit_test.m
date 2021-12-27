close all
 tic
format % long
clearvars -except  yt  peter_Swain_data
%fix_noise=false;%true;
compare_to_prevus=0; % 0 not work !
%figure(10)
hold on
kkkk=zeros(19,400);
namber_of_experiments=6;

%for j=1:namber_of_experiments %123
    for j=14:19
j


f=1;
nf=1;



yt=peter_Swain_data(j).yt;
%yt=round(yt,2);
%yt_big=yt;
%yt = zip_yt( yt,5 );



f_name=peter_Swain_data(j).name;
% start_frame=6; %123
% end_frame=36;%40;
start_frame=2;
end_frame=20;%40;

frame_num=end_frame-start_frame+1;
[cell_num,~]=size(yt);
x=zeros(cell_num,frame_num);
xx=zeros(cell_num,frame_num);
for i=1:cell_num
    
    
    
    y=yt(i,start_frame:end_frame)';
   %y=yt_big(i,:)';
    %[fitresult, gof] = two_exp_fit(y);
    [fitresult, gof] = two_exp_fit(y);
    %y=yt(i,start_frame:end_frame)';
   % if gof.rsquare>0.96
   kkkk(j,i)=gof.rsquare;
   if gof.rsquare>0.95
        %ey=feval(fitresult, 1:frame_num);
        ey=feval(fitresult, 1:frame_num);
        
        if compare_to_prevus==1
            p=ey(2:end)./ey(1:end-1);
            ey(2:end)=y(1:end-1).*p;
            
         end   
            
            
            
       
    
    
    y_noise_fix=zeros(1,length(y));
    
    for k=1:length(y)
        if (y(k)-ey(k))^2>ey(k)
        y_noise_fix(k)=(y(k)-ey(k))^2-ey(k);
         f=f+1;
         else
            y_noise_fix(k)=(y(k)-ey(k))^2;
            nf=nf+1;
         end
        
    end
   y_noise_fix=y_noise_fix';
    y_noise=(y-ey).^2;

    
       
    % y_noise_fix=y_noise_fix';
    
    
    
    
    
    
    
    
    x(i,1:frame_num)=ey.*(y(1)-ey)./y_noise;
   xx(i,1:frame_num)=ey.*(y(1)-ey)./y_noise_fix;
    
    
    end
end

% f
% nf

% x(x<=0)=[];
% xx(xx<=0)=[];
x=real(x);
xx=real(xx);

v=log10(real(x(:)));
v(v<-10000000000000)=[];
v=real(v);

%subplot(1,namber_of_experiments,j) %123
subplot(1,namber_of_experiments,j-13)
hold on
histogram(real(v))
mean(xx(:))
vv=log10(real(xx(:)));
vv(vv<-10000000000000)=[];

vv=real(vv);

histogram(real(vv))
legend('standart','fix')
title([f_name ' ' num2str(median(v(:))) ' ' num2str(median(vv(:)))])
xlim([1 8])

res(j,1)=median(v(:));
res(j,2)=median(vv(:));
res(j,3)=var(v(:));
res(j,4)=var(vv(:));

end


for d=1:namber_of_experiments
    %t(d)=peter_Swain_data(d).EX_19_res;
    t(d)=peter_Swain_data(d).true_res;
end
toc
res(:,1)./t'
res(:,2)./t'
 toc
 
 res(:,1)-t'
res(:,2)-t'