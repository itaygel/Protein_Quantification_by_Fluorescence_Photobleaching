format long
close all
figure(10)
hold on

clearvars -except  peter_Swain_data all_data
ex_num=6
c=[]
for j=1:ex_num

yt=peter_Swain_data(j).yt;
f_name=peter_Swain_data(j).name;
ey=peter_Swain_data(j).ey;

[a, b]=size(ey);
% for e=1:a
%    ey(e,1:b-1)=smooth(ey(e,1:b-1),10); 
%  
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  peter_Swain_data(13+j).name=['f_name',' zip 2'];
%  peter_Swain_data(13+j).true_res=peter_Swain_data(j).true_res;
%  peter_Swain_data(13+j).yt=zip_yt( yt,2 );
% 
% c=[c;peter_Swain_data(13+j).yt];

w=all_data(1:length(yt),:);
all_data(1:length(yt),:)=[];
peter_Swain_data(13+j).ey=w;
%%%%%%%%%%%%%%%%%%%%%%%%%%



l=length(yt);
if j<=7
    d=45;
else
    if j<=13
    d=11;
    else
        d-22;
    end
end
for i=3:10
%     eyy=ey(:,i);
%    x(i-1,1:l)=eyy.*(yt(:,1)-eyy)./((yt(:,i)-eyy).^2);
    %x(i-1,1:l)=ey(:,i).*(yt(:,1)-ey(:,i))./((yt(:,i)-ey(:,i)).^2);
%    plot(1:l, yt(:,i),'.r')
%    plot( (1:l)+0.2,ey(:,i),'.g')yt=1
 %xx(i-1,1:l)=ey(:,i).*(yt(:,1)-ey(:,i))./((yt(:,i)-ey(:,i)).^2-yt(:,i));
p=ey(:,i)./ey(:,1);
pp=p.*(1-p);
ni(i-2,1:l)=((yt(:,i)-ey(:,i)).^2)./(yt(:,1).*pp);
x(i-2,1:l)=yt(:,1)./ni(i-2,1:l)';

ni2(i-2,1:l)=((yt(:,i)-ey(:,i)).^2-yt(:,1))./(yt(:,1).*pp);
xx(i-2,1:l)=yt(:,1)./ni2(i-2,1:l)';

p=ey(:,i)./ey(:,i-1);
pp=p.*(1-p);
%ni3(i-2,1:l)=((yt(:,i)-ey(:,i)).^2)./(yt(:,i-1).*pp);
ni3(i-2,1:l)=((yt(:,i)-yt(:,i-1).*p).^2)./(yt(:,i-1).*pp);
xxx(i-2,1:l)=yt(:,1)./ni3(i-2,1:l)';


end
% 

xx=xxx;

subplot(1,ex_num,j)
hold on
% xx=x(:);
% xx(xx>10^5)=[];
% xx(xx<0)=[];
% mean(xx)
v=log10(x(:));
vv=log10(xx(:));
histogram(real(v))
hold on
histogram(real(vv))
legend('standart','fix')

%histogram(xx,100,'BinLimits',[20,10^4])
legend('standart','fix')
title([f_name ' ' num2str(median(v(:))) ' ' num2str(median(vv(:)))])
xlim([1 8])

res(j,1)=median(v(:));
res(j,2)=median(vv(:));
res(j,3)=var(v(:));
res(j,4)=var(vv(:));
end



for d=1:ex_num
    t(d)=peter_Swain_data(d).true_res;
end
res(:,1)./t'
res(:,2)./t'