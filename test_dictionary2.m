%% 1.18  测试KSVD和JBHI过完备字典
% 对比KSVD和JBHI方法的误差和耗时 
%19种压缩比 ，对100段ECG取平均值
% 重构算法固定bp，测量矩阵固定一个伯努利矩阵
% JBHI方法的QRS检测方法经过改进
%%                                
load D_S.mat    %基本字典
load D_NQp.mat  %不含QRS的子字典
load D_Qp.mat      %论文使用子字典（12个）
N=192; 
q=12;
% threshold=1.5;     %判断QRS的阈值
load('BernoulliSample.mat');   %导入一个固定测量矩阵 1024*1024

for num_CR=1:19          %19种CR
    fprintf('>>压缩比%d，开始测试\n',num_CR)
    M=fix((1-num_CR*5*0.01)*N);             %计算CR对应测量数  CR从95%、90% ~ 5%，对应M=972、921~51   
    Phi=BernoulliSample(1:M,1:N);                %每种CR对应一种测量数M，对应一种测量矩阵M*N     
      
                  
                  fprintf('>压缩比%d，开始测试\n',num_CR)
                  toc_sum=0;      %计时器和归零
                  PRDsum_KSVD=0;       %误差和清零  （为了计算一种算法在一种CR下  对100段心电下误差的平均值
                  PRDsum_JBHI=0; 
                  TIMEsum_KSVD=0;
                  TIMEsum_JBHI=0;
                  
                   for num_ECG=1:100;   
                            fprintf('loading ECG%d\n',num_ECG)
                            ecgstr=['ecg',num2str(num_ECG)];                 
                            x=cell2mat(struct2cell(load(ecgstr)));        %每次导入一个ecg  从ecg1~ecg100
                            x_dic=x(1:960);
                         
                            %%   KSVD方法
                            L=5;       %信号数量
                            A=reshape(x_dic,192,5);
                            X=[ ];                      %重构后 将各小段连接的一维信号 
                            tic;
                            for i=1:L
                                x=A(:,i); 
                                %Phi=BernoulliMtx( M,N );                                          
                                y=Phi*x;                                                             
                                Psi=D_S;
                                T=Phi*Psi; 
                                % hat_s=cs_irls(y,T,N);                     %迭代重加权
                                hat_s=cs_bp(y,T,N);                         %基追踪
                                hat_x=real(Psi*hat_s);
                                X=[X hat_x']; 
                            end
                            X_s=X';                                           %X_s是使用标准字典的重构信号
                            TIME_KSVD=toc;
                            PRD_KSVD=norm(x_dic-X_s)/norm(x_dic)*100;     % 标准字典的重构误差
                            
                           %% JBHI方法        
region_n=12;            %子区域数量
region_l=N/region_n; %区域宽度           
win_w=N/q;    %窗口宽度
Q=cell(1,q);             %元胞数组q个
X_sub=[];
num=5;
tic;
[R,V,R_no]=QRS_crest(x_dic); 
R_num=length(R);
for i=1:num
%% 对信号进行QRS区域判定
    p_start=N*(i-1)+1;
    p_end=N*i;
    x=x_dic(p_start:p_end);
    flag_qrs=0;
        for ii=1:R_num
            if R(ii)>=p_start&&R(ii)<=p_end
                ww=ceil((R(ii)-p_start+1)/region_l);  %w为波峰所在窗口区域 位置在1~q中的一个
                dictionary=D_Qp{ww};   %使用含QRS的字典
                flag_qrs=1;    
            end     
        end
        
        if  flag_qrs==0
            dictionary=D_NQp;    %使用不含QRS的字典
        end
        
% Phi=BernoulliMtx( M,N );                                          
y=Phi*x;                                                             
Psi=dictionary;
T=Phi*Psi; 
% hat_s=cs_irls(y,T,N);                     %迭代重加权
hat_s=cs_bp(y,T,N);                         %基追踪
hat_x=real(Psi*hat_s);
X_sub=[X_sub hat_x'];                   %分段信号拼接
end
X_sub=X_sub';
TIME_JBHI=toc;
TIME_JBHI=TIME_JBHI+TIME_KSVD;
PRD_JBHI=norm(x_dic-X_sub)/norm(x_dic)*100;  % 子字典的重构误差  
%% 统计总误差和总时间
PRDsum_KSVD=PRDsum_KSVD+PRD_KSVD;
PRDsum_JBHI=PRDsum_JBHI+PRD_JBHI;
TIMEsum_KSVD=TIMEsum_KSVD+TIME_KSVD;
TIMEsum_JBHI=TIMEsum_JBHI+TIME_JBHI;

                    end           %  100段ECG测试完成
                    
                    PRDaver_KSVD=PRDsum_KSVD/num_ECG;
                    PRDaver_JBHI=PRDsum_JBHI/num_ECG;             
                    TIMEaver_KSVD=TIMEsum_KSVD/num_ECG;
                    TIMEaver_JBHI= TIMEsum_JBHI/num_ECG;
                    
                    Sheet1(1,num_CR)=PRDaver_KSVD;     % KSVD 误差
                    Sheet1(2,num_CR)=PRDaver_JBHI;      % JBHI  误差
                    Sheet2(1,num_CR)=TIMEaver_KSVD;     %KSVD 耗时     
                    Sheet2(2,num_CR)=TIMEaver_JBHI;     %JBHI  耗时

                    fprintf('>压缩比%d，完成测试 \n',num_CR)
                   
          fprintf('>>压缩比%d，完成所有测试 \n',num_CR)
end          % 19种CR 遍历结束
fprintf('>>>完成所有测试，正将数据存入表格... \n')

xlswrite('KSVD&JBHI_19CR_100ecg.xlsx',Sheet1,'Sheet1','B3');
xlswrite('KSVD&JBHI_19CR_100ecg.xlsx',Sheet2,'Sheet2','B3');
fprintf('>>>>>>>> 已保存数据到表格<<<<<<<<\n')
sprintf('>>All Completed<<\n')


% figure(1);
% hold on;
% plot(Sheet1(1,:),'r')
% plot(Sheet1(2,:),'g') 
% legend('KSVD','JBHI')  
% 
% 
% figure(2);
% plot(Sheet2(1,:),'r')
% plot(Sheet2(2,:),'g') 
% legend('KSVD','JBHI')  





