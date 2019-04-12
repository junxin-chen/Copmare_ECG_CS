%************************************************************************%
function hat_x=cs_sp(y,T_Mat,m)

% Reference: W. Dai and O. Milenkovic, ¡°Subspace Pursuit for Compressive
% Sensing Signal Reconstruction,¡± 2009.


% y=T_Mat*x, T_Mat is n-by-m
% y - measurements
% T_Mat - combination of random matrix and sparse representation basis
% m - size of the original signal
% the sparsity is length(y)/4

n=length(y);                           % length of measurements
s=floor(n/4);                                 % sparsity                  
r_n=y;                                 % initial residuals

sig_pos_lt=[];                         % significant pos for last time iteration

for times=1:s                          % number of iterations
    
    product=abs(T_Mat'*r_n);
    [val,pos]=sort(product,'descend');
    sig_pos_cr=pos(1:s);             % significant pos for curretn iteration
    
    sig_pos=union(sig_pos_cr,sig_pos_lt);
    
    Aug_t=T_Mat(:,sig_pos);            % current selected entries of T_Mat 
    
    aug_x_cr=zeros(m,1);               
    aug_x_cr(sig_pos)=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;  % temp recovered x (sparse)
    
    [val,pos]=sort(abs(aug_x_cr),'descend');
    
    hat_x=zeros(1,m);
    hat_x(pos(1:s))=aug_x_cr(pos(1:s));% recovered x with s sparsity  
    
    sig_pos_lt=pos(1:s);               % refresh the significant positions
    
    r_n=y-T_Mat*hat_x';
end
             














