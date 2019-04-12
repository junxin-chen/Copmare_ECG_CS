%% Tompkins差分法检测QRS波，找到R波波峰，返回每个波峰的横坐标

function [R,V,R_no]=QRS_crest(x)
%R为波峰横坐标；V为波峰纵坐标；R_no为QRS波段数

m=length(x);   % =SAMPLES2READ
y1 = [];%一阶差分
for n = 1:(m-1)
    y1(n) = x(n+1) - x(n);
end
y1(m)=y1(m-1);

y2 = [];%二阶差分
for n = 2:(m-1)
    y2(n) = x(n+1)-2*x(n)+x(n-1);
end
y2(1)=y2(2);
y2(m)=y2(m-1);

y=[]; %Tompkins差分法是取心电信号的一阶与二阶差分的平方和作为QRS波群输出标记的脉冲信号
for k=1:1:m;
    y(k)=y1(k)^2+y2(k)^2;
end

 th=0.02;  % 经验值,可回溯检查,调整阈值--处理后的R波幅值下限
 QRS=find(y-th>0);
 qrs={};  % QRS波群
 
 r=1;
 j=2;
 qrs{1}(1)=QRS(1);%qrs{r}用于记录R波时间坐标
 
for i=2:1:length(QRS);
    if QRS(i)-QRS(i-1)<=100;%时间坐标相差100以下认为没有到下个QRS波段，将这个QRS段时间坐标存入qrs{r}
        qrs{r}(j)=QRS(i);
            j=j+1;
    else
        qrs{r+1}(1)=QRS(i);%时间坐标相差100以上认为到下个QRS波段，将这个QRS段时间坐标存入qrs{r+1}
        r=r+1;
        R_no=r;%记录QRS波段数
        j=2;   
    end 
end

Rsite=zeros(R_no,1);%QRS段数为行数，10列，每个QRS波段记录10个参数


for r=1:1:R_no;
    n=1;
    for i=2:1:length(qrs{r});
    if qrs{r}(i)-qrs{r}(i-1)>1;%如果此记录点的时间坐标比前一个大1以上
        Rsite(r)=qrs{r}(i-1)+2;%将前一个记录点的时间坐标加2记录到Rsite
        break;
    end
    end
    if Rsite(r)==0;
        Rsite(r)=min(qrs{r});
    end
end

R=zeros(R_no,1);%R波波锋的横坐标
V=zeros(R_no,1);%R波波锋的纵坐标

for i=1:R_no;
[V(i),R(i)]=max(x(Rsite(i)-20:Rsite(i)+20));
R(i)=R(i)+Rsite(i)-21;
end
