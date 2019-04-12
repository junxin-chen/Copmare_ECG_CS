%************************************************************************%
function hat_x=cs_omp(y,T_Mat,m)
% Reference: J. Tropp and A. Gilbert, “Signal Recovery from Random 
% Measurements via Orthogonal Matching Pursuit,” 2007.

% y=T_Mat*x, T_Mat is n-by-m
% y - measurements
% T_Mat - combination of random matrix and sparse representation basis
% m - size of the original signal
% the sparsity is length(y)/4

n=length(y);
s=4*floor(n/4);                                     %  测量值维数
hat_x=zeros(1,m);                                 %  待重构的谱域(变换域)向量                     
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_n=y;                                            %  残差值 

for times=1:s;                                  %  迭代次数(稀疏度是测量的1/4)

    product=abs(T_Mat'*r_n);
    
    [val,pos]=max(product);                       %  最大投影系数对应的位置
    Aug_t=[Aug_t,T_Mat(:,pos)];                   %  矩阵扩充
    T_Mat(:,pos)=zeros(n,1);                      %  选中的列置零（实质上应该去掉，为了简单将其置零）
    aug_x=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;           %  最小二乘,使残差最小
    r_n=y-Aug_t*aug_x;                            %  残差
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
    
end
hat_x(pos_array)=aug_x;                           %  重构的向量 


