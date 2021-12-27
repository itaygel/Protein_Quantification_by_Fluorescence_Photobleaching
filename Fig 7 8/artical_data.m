format long
close all
figure(10)
hold on

clearvars -except  peter_Swain_data dd dd1
ex_num=18;


all_sn_noise=[];
all_pb_noise=[];



for j=1:ex_num
x=[];
xx=[];
y_noise_fix=[];
ey_ctp=[];

yt=peter_Swain_data(j).yt;
f_name=peter_Swain_data(j).name;
ey=peter_Swain_data(j).ey;
fit_x=peter_Swain_data(j).ex19_plot(:,1);
fit_y=peter_Swain_data(j).ex19_plot(:,2);
[a, b]=size(ey);
% for e=1:a
%    ey(e,1:b-1)=smooth(ey(e,1:b-1),10); 
%  
% end
%%%%%%%%%%%%%%%%%%%%%%%%%
% p=ey(:,[2:end])./ey(:,[1:end-1]);
% y0=yt(:,1);
%             ey_ctp=ey;      
%             ey_ctp([1:a],[2:end])=yt(:,[1:end-1]).*p; 
%              y_noise_ctp=(yt-ey_ctp).^2;
%            % xxx([1:a-1],1:frame_num-1)=(y0./yt(2:end)).*ey_ctp(2:end).*(yt(1:end-1)-ey_ctp(2:end))./y_noise_ctp(2:end);
% xxx_pfc=(y0./yt(2:end)).*ey_ctp(2:end).*(yt(1:end-1)-ey_ctp(2:end))./y_noise_ctp(2:end);
           %%%%%%%%%%%%%%%%%%%%%%%%%

if j<=6
    d=45;
   % d=40;
else
    if j<=12
    d=22;
    %d=20;
    else
        d=11;
        %d=10;
    end
end

%d=45;
l=length(yt);
for i=2:d
%     eyy=ey(:,i);
%    x(i-1,1:l)=eyy.*(yt(:,1)-eyy)./((yt(:,i)-eyy).^2);
    x(i-1,1:l)=ey(:,i).*(yt(:,1)-ey(:,i))./((yt(:,i)-ey(:,i)).^2);
   %x2(i-1,1:l)=ey(:,i).*(yt(:,1)-ey(:,i))./((yt(:,i)-ey(:,i)).^2);

%    plot(1:l, yt(:,i),'.r')
%    plot( (1:l)+0.2,ey(:,i),'.g')yt=1
 %xx(i-1,1:l)=ey(:,i).*(yt(:,1)-ey(:,i))./((yt(:,i)-ey(:,i)).^2-yt(:,i));
 
  for k=1:length(yt)
        %if (yt(k,i)-ey(k,i))^2>ey(k,i)
        if 2>1
        y_noise_fix(k)=(yt(k,i)-ey(k,i))^2-ey(k,i);
        % f=f+1;
         else
            y_noise_fix(k)=(yt(k,i)-ey(k,i))^2;
        %    nf=nf+1;
         end
        
  end
    
 % y_noise_fix=(yt-ey).^2-ey;
 

 xx(i-1,1:l)=ey(:,i).*(yt(:,1)-ey(:,i))./y_noise_fix';
 
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%
%  all_sn_noise=[all_sn_noise;yt(:,i)];
%  all_pb_noise=[all_pb_noise;(yt(:,i)-ey(:,i)).^2];
 %%%%%%%%%%%%%%%%%%%%%%%
end
% 


%figure(j)
subplot(1,ex_num,j)
xlim([1 8])
hold on
% xx=x(:);
% xx(xx>10^5)=[];
% xx(xx<0)=[];
% mean(xx)
v=log10(x(:));
vv=log10(xx(:));


%vvv_pfc=log10(xxx_pfc(:));
%[a]=histogram(real(v),'Normalization','probability','BinWidth',1);
histogram(real(v),'LineWidth',0.5,'Normalization','probability');

hold on
%histogram(real(vv),'Normalization','probability','BinWidth',1)

%histogram(real(vv),'LineWidth',0.5,'Normalization','probability');
%plotyy(nan,nan,fit_x,fit_y)

%legend('standart','fix')

%plot_19ex_data(peter_Swain_data(j).ex19_plot,a.Values)

%histogram(xx,100,'BinLimits',[20,10^4])
%legend('standart','fix')
%title([f_name ' ' num2str(median(v(:))) ' ' num2str(median(vv(:)))])
title([f_name])


res(j,1)=median(v(:));
res(j,2)=median(vv(:));
res(j,3)=var(v(:));
res(j,4)=var(vv(:));
%res(j,5)=var(vvv_pfc(:));
%res(j,6)=var(vvv_pfc(:));


end



for d=1:ex_num
    t(d)=peter_Swain_data(d).true_res;
end

res(:,2)./t'

a=res(:,1)./t';





figure(1001)

hold on
plot(t(1:6),res(1:6,1),'ob')
plot(t(1:6),res(1:6,2),'xr')
plot([3 6.5],[3 6.5],'k:')
errorbar(t(1:6),res(1:6,1),sqrt(res(1:6,3)),'b','LineStyle','none')
errorbar(t(1:6),res(1:6,2),sqrt(res(1:6,4)),'r','LineStyle','none')

legend('no fix','fix')
%legend('FFC','PFC')

xlabel('Biochemical experiments  log_{10}(number of molecules)')
ylabel('Bi-exponential fit log_{10}(number of molecules)')















figure(234)

plot(a(1:6,1),a(7:12,1),'x')
hold on
plot(a(1:6,1),a(13:18,1),'ob')
legend('8-13','14-19')







figure(2354)
%plot(res(1:6,1)./res(14:19,1),'xb')
hold on
%plot(res(1:6,1)./res(8:13,1),'xk')
plot(t(1:6)'./res(1:6,1),'*g')
plot(t(1:6)'./res(7:12,1),'*r')
plot(t(1:6)'./res(13:18,1),'*b')

legend('1-6','8-13','14-19')



figure(5)  
%plot(res(1:6,1)./res(14:19,1),'xb')
hold on


plot(t(1:6),res(1:6,1),'or')
%plot(t(1:6),res(1:6,2),'xr')
%plot([3 6.5],[3 6.5],'k')
%errorbar(t(1:6),res(1:6,2),sqrt(res(1:6,3)),'b','LineStyle','none')
%errorbar(t(1:6),res(1:6,7),sqrt(res(1:6,6)),'r','LineStyle','none')


xlabel('Biochemical experiments')
ylabel('Estimates from bi-exponential fit')




plot(t(1:6),res(7:12,1),'xb')
%plot(t(1:6),res(7:12,2),'xb')

plot(t(1:6),res(13:18,1),'+k')




figure(100)
for i=1:2
    subplot(2,1,i) 
    a=[1,2,2,7];
    b=a(i);
%plot(res(1:6,1)./res(14:19,1),'xb')
hold on
%plot(res(1:6,1)./res(8:13,1),'xk')
%plot(t(1:6)'./res(1:6,1),'*g')

%plot(res(7:12,b)./res(1:6,b),'*r')
%plot(res(13:18,b)./res(1:6,b),'*b')

plot(10.^(res(7:12,b))./10.^(res(1:6,b)),'*r')
plot(10.^(res(13:18,b))./10.^(res(1:6,b)),'*b')
legend('8-13','14-19')

end



figure(777)
%subplot(4,1,1)
e=1;
plot(t(1:6),res([1:6],e),'*')
hold on
plot(t(1:6),res([7:12],e),'*')
plot(t(1:6),res([13:18],e),'*')
%plot(t(1:6),t(1:6),'+')
plot([3.5 6],[3.5 6],':')

legend('\mu','2\mu','4\mu')
xlabel('Biochemical experiments')
ylabel('Gaussian Processes Estimating')


figure(778)
%subplot(4,1,1)
e=1;
plot(10.^(t(1:6)),10.^(res([1:6],e)),'*')
hold on
plot(10.^(t(1:6)),10.^(res([7:12],e)),'*')
plot(10.^(t(1:6)),10.^(res([13:18],e)),'*')
%plot(t(1:6),t(1:6),'+')
plot(10.^([3.5 6]),10.^([3.5 6]),':')

legend('\mu','2\mu','4\mu')
xlabel('Biochemical experiments')
ylabel('Gaussian Processes Estimating')

ey=peter_Swain_data(1).ey;
ey=ey(:,2:end-1)./ey(:,1:end-2);
figure()
%plot(ey(20,:))
hold on
%plot(ey(26,:))
%plot(ey(79,:))

plot(ey(1,:))
plot(ey(2,:))
plot(ey(3,:))
plot(ey(4,:))
plot(ey(5,:))
%
% for i=1:4
% figure(7737)
% %subplot(4,1,1)
% a=[1 2 3 4 9];
% e=a(i);
% subplot(1,4,i)
% plot(t(1:6),res([1:6],e),'o')
% hold on
% plot(t(1:6),res([7:12],e),'o')
% plot(t(1:6),res([13:18],e),'o')
% %plot(t(1:6),t(1:6),'+')
% plot([3.5 6],[3.5 6],':')
% legend('\mu','2\mu','4\mu')
% xlabel('Biochemical experiments')
% ylabel('Gaussian Processes Estimating')
% end


%set(gcf,'PaperUnits','cent','PaperPosition',[0 0 18.5 12])
%print(['6subplot'], '-dtiff', '-r900')
