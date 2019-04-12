%% 1.18  测试小波基
% 对比52种小波基的误差 和耗时
% 重构算法固定bp，测量矩阵固定伯努利矩阵
% 19种压缩比 ，对100段ECG取平均值

%%                                
N=1024;   
load('BernoulliSample.mat');   %导入一个固定测量矩阵 1024*1024
for num_CR=1:19          %19种CR
    fprintf('>>压缩比%d，开始测试\n',num_CR)
    M=fix((1-num_CR*5*0.01)*N);             %计算CR对应测量数  CR从95%、90% ~ 5%，对应M=972、921~51   
    Phi=BernoulliSample(1:M,:);                %每种CR对应一种测量数M，对应一种测量矩阵M*N
   % Phi=BernoulliMtx( M,N );        
         for num_wave=1:52        % 52种小波基 
                   fprintf('>压缩比%d-小波基%d，开始测试\n',num_CR,num_wave)
                   toc_sum=0;      %计时器和归零
                   PRD_sum=0;       %误差和清零  （为了计算一种算法在一种CR下  对100段心电下误差的平均值
                    
                   for num_ECG=1:100;    %100段ECG
                            fprintf('loading ECG%d\n',num_ECG)
                            ecgstr=['ecg',num2str(num_ECG)];                 
                            x=cell2mat(struct2cell(load(ecgstr)));        %每次导入一个ecg  从ecg1~ecg100
                            y=Phi*x; 
                            tic       %开始计时   重建信号耗时
                            [ww]=dwtmtxx( N, num_wave,5);    %小波基的选择  在函数dwtmtxx.m中
                            Psi=[ww];                                           
                            T=Phi*Psi'; 
                            hat_s=cs_bp(y,T,N);                     %基追踪算法BP
                            hat_x=real(Psi'*hat_s);
                            time_end=toc;   %结束计时   重建信号耗时
                            toc_sum=toc_sum+time_end;          %累计时间
                            PRD=norm(x-hat_x)/norm(x)*100;     %计算本次误差
                            PRD_sum=PRD_sum+PRD;                %累计误差      
                    end           %  100段ECG测试完成
                    
                    fprintf('>压缩比%d-小波基%d，完成测试 \n',num_CR,num_wave)
                    PRD_aver=PRD_sum/num_ECG;  %100段ECG的误差平均值
                    time_aver=toc_sum/num_ECG;   %100段ECG的耗时平均值
                    A(num_wave,num_CR)=PRD_aver;     %A表 误差数据
                    B(num_wave,num_CR)=toc_sum;     %B表 耗时数据         
          end        %一种压缩比对应的小波基遍历结束
    
          fprintf('>>压缩比%d，完成所有测试 \n',num_CR)
end          % 19种CR 遍历结束
fprintf('>>>完成所有测试，正将数据存入表格... \n')

xlswrite('sparse_19CR_100ecg.xlsx',A,'Sheet1','B3');
xlswrite('sparse_19CR_100ecg.xlsx',B,'Sheet2','B3');
fprintf('>>>>>>>> 已保存数据到表格<<<<<<<<\n')
sprintf('>>All Completed<<\n')







