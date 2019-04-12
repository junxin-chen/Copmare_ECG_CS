%% 1.18  测试重构算法。
%对比5种算法的误差和耗时
%稀疏基固定小波基db2，测量矩阵固定伯努利矩阵 
%19种压缩比  对100段ECG结果取平均值

%%                                
N=1024;   
load('BernoulliSample.mat');   %导入一个固定测量矩阵 1024*1024
for num_CR=1:19          %19种CR
    fprintf('>>压缩比%d，开始测试\n',num_CR)
    M=fix((1-num_CR*5*0.01)*N);             %计算CR对应测量数  CR从95%、90% ~ 5%，对应M=972、921~51   
    Phi=BernoulliSample(1:M,:);                %每种CR对应一种测量数M，对应一种测量矩阵M*N
   % Phi=BernoulliMtx( M,N );        
         for num_algo=1:5       % 5种重构算法 
                   fprintf('>压缩比%d-算法%d，开始测试\n',num_CR,num_algo)
                   toc_sum=0;      %计时器和归零
                   PRD_sum=0;       %误差和清零  （为了计算一种算法在一种CR下  对100段心电下误差的平均值
                    for num_ECG=1:100;    %100段ECG
                            fprintf('loading ECG%d\n',num_ECG)
                            ecgstr=['ecg',num2str(num_ECG)];                 
                            x=cell2mat(struct2cell(load(ecgstr)));        %每次导入一个ecg  从ecg1~ecg100
                            y=Phi*x;                                                     
                            [ww]=dwtmtx( N,'db2',5);     %选择小波基
                            Psi=[ww];                                  
                            T=Phi*Psi'; 
                            tic       %开始计时   重建信号耗时
                                        switch num_algo                                         %重构算法选择  K决定
                                            case 1
                                                 hat_s=cs_omp(y,T,N);
                                                 hat_x=real(Psi'*hat_s.');                  
                                            case 2
                                                hat_s=cs_bp(y,T,N);    
                                                hat_x=real(Psi'*hat_s);        
                                            case 3
                                                 hat_s=cs_cosamp(y,T,N);
                                                 hat_x=real(Psi'*hat_s.');         
                                            case 4
                                                 hat_s=cs_irls(y,T,N);  
                                                 hat_x=real(Psi'*hat_s);       
                                            case 5
                                                 hat_s=cs_sp(y,T,N); 
                                                 hat_x=real(Psi'*hat_s.');       
                                        end
                                time_end=toc;   %结束计时   重建信号耗时
                                toc_sum=toc_sum+time_end;          %累计时间
                                PRD=norm(x-hat_x)/norm(x)*100;     %计算本次误差
                                PRD_sum=PRD_sum+PRD;                %累计误差      
                                                                   
                    end           %  100段ECG测试完成
                    fprintf('>压缩比%d-算法%d，完成测试 \n',num_CR,num_algo)
                    PRD_aver=PRD_sum/num_ECG;  %100段ECG的误差平均值
                    time_aver=toc_sum/num_ECG;   %100段ECG的耗时平均值

                    A(num_algo,num_CR)=PRD_aver;     %A表 误差数据
                    B(num_algo,num_CR)=toc_sum;     %B表 耗时数据
                   
          end        %重构算法遍历结束
          fprintf('>>压缩比%d，完成测试 \n',num_CR)
end          % 19种CR 遍历结束
fprintf('>>>完成所有测试，正将数据存入表格... \n')

xlswrite('algo_19CR_100ecg.xlsx',A,'Sheet1','B3');  %生成重构算法误差对比的表格
xlswrite('algo_19CR_100ecg.xlsx',B,'Sheet2','B3');  %生成重构算法耗时对比的表格
fprintf('>>>>>>>> 已保存数据到表格<<<<<<<<\n')
sprintf('>>All Completed<<\n')







