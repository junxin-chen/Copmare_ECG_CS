function [ ww ] = dwtmtx( N,w,wlev )  
%DWTMTX Discrete wavelet transform matrix  
%   This function generates the transform matrix ww according to input   此函数根据输入生成变换矩阵WW
%   parameters N,wtype,wlev .  
%Detailed explanation goes here  
%   N is the dimension of ww    ww的维度
%   wtype is the wavelet type  小波类型
%   wlev is the number of decomposition level  分解的层数
%NOTE: The extension mode must be Periodization('per')  



wav(1)={'haar'};
wav(2)={'db2'};
wav(3)={'db3'};
wav(4)={'db4'};
wav(5)={'db5'};
wav(6)={'db6'};
wav(7)={'db7'};
wav(8)={'db8'};
wav(9)={'db9'};
wav(10)={'db10'};
wav(11)={'sym2'};
wav(12)={'sym3'};
wav(13)={'sym4'};
wav(14)={'sym5'};
wav(15)={'sym6'};
wav(16)={'sym7'};
wav(17)={'sym8'};
wav(18)={'coif1'};
wav(19)={'coif2'};
wav(20)={'coif3'};
wav(21)={'coif4'};
wav(22)={'coif5'};
wav(23)={'bior1.1'};
wav(24)={'bior1.3'};
wav(25)={'bior1.5'};
wav(26)={'bior2.2'};
wav(27)={'bior2.4'};
wav(28)={'bior2.6'};
wav(29)={'bior2.8'};
wav(30)={'bior3.1'};
wav(31)={'bior3.3'};
wav(32)={'bior3.5'};
wav(33)={'bior3.7'};
wav(34)={'bior3.9'};
wav(35)={'bior4.4'};
wav(36)={'bior5.5'};
wav(37)={'bior6.8'};
wav(38)={'rbio1.1'};
wav(39)={'rbio1.3'};
wav(40)={'rbio1.5'};
wav(41)={'rbio2.2'};
wav(42)={'rbio2.4'};
wav(43)={'rbio2.6'};
wav(44)={'rbio2.8'};
wav(45)={'rbio3.1'};
wav(46)={'rbio3.3'};
wav(47)={'rbio3.5'};
wav(48)={'rbio3.7'};
wav(49)={'rbio3.9'};
wav(50)={'rbio4.4'};
wav(51)={'rbio5.5'};
wav(52)={'rbio6.8'};


wtype=wav{w};

[h,g]= wfilters(wtype,'d');         %Decomposition low&high pass filter  分解低通高通滤波器
L=length(h);                        %Filter length  滤波器长度
h_1 = fliplr(h);                    %Flip matrix left to right  从左到右翻转矩阵
g_1 = fliplr(g);  
loop_max = log2(N);  
loop_min = double(int8(log2(L)))+1;  
if wlev>loop_max-loop_min+1  
    fprintf('\nWaring: wlev is too big\n');  
    fprintf('The biggest wlev is %d\n',loop_max-loop_min+1);  
    wlev = loop_max-loop_min+1;  
end  
ww=1;  
for loop = loop_max-wlev+1:loop_max  
    Nii = 2^loop;  
    p1_0 = [h_1 zeros(1,Nii-L)];  
    p2_0 = [g_1 zeros(1,Nii-L)];  
    p1 = zeros(Nii/2,Nii);  
    p2 = zeros(Nii/2,Nii);  
    for ii=1:Nii/2  
        p1(ii,:)=circshift(p1_0',2*(ii-1)+1-(L-1)+L/2-1)';  
        p2(ii,:)=circshift(p2_0',2*(ii-1)+1-(L-1)+L/2-1)';  
    end  
    w1=[p1;p2];  
    mm=2^loop_max-length(w1);  
    w=[w1,zeros(length(w1),mm);zeros(mm,length(w1)),eye(mm,mm)];  
    ww=ww*w;  
    clear p1;clear p2;  
end  
  
%The end!!!  
end  