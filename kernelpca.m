function [eig_vector, lambdas] = kernelpca(data_in,var,type,par1,par2)
%{ 
eig_vector (OUTPUT) = First k Eigenvectors (Principal Components)    
lambdas (OUTPUT)=  First k eigenvalues 
var = % of Variance we want to retain in the transformed data
type = Choose the kernel function to use: 'G' for Gaussian, 'P' for Polynomial, 'S' for Sigmoid
par1, par2 = values of parameters for the chosen kernel (look at the definition of each type of kernel function
             below to check what is being represented by par1 and par2)   
Data_in - Input data (d (dimensions) X N (# of points) 
%}
N = normalize(data_in');
data_in = N';
%% STEP 1
%% Using the Gaussian Kernel to construct the Kernel K
% K(x,y) = -exp(par1*||x-y||^2) 
K = zeros(size(data_in,2),size(data_in,2));
if type == 'G'
     sq_dists = pdist(data_in', 'squaredeuclidean');
     mat_sq_dists = squareform(sq_dists);
     K = exp(-par1 * mat_sq_dists);    
end
%% Using the Polynomial Kernel to construct the Kernel K
% K(x,y) = (xi, xj + k)^d
if type == 'P'
for row = 1:size(data_in,2)
    for col = 1:row
      temp = (dot(data_in(:,row), data_in(:,col))+ par2)^par1;
      K(row,col) = temp;
    end
end
K = K + K'; 
% Dividing the diagonal element by 2 since it has been added to itself
for row = 1:size(data_in,2)
    K(row,row) = K(row,row)/2;
end
end
%% Using the Sigmoid Kernel to construct the Kernel K
% K(x,y) = tanh(alpha*(xT.y) + c)
if type == 'S'
for row = 1:size(data_in,2)
    for col = 1:row
      temp = tanh(par1*(data_in(:,row)'*data_in(:,col)) + par2);
      K(row,col) = temp;
    end
end
K = K + K'; 
% Dividing the diagonal element by 2 since it has been added to itself
for row = 1:size(data_in,2)
    K(row,row) = K(row,row)/2;
end
end


%% STEP 2
%{
We know that for Kernel PCA the data has to be centered. Even if the input data set is centered, 
there is no gurantee that the data when mapped in the feature space is also centered. 
Since we actually never work in the feature space we cannot center the data. 
To include this correction a pseudo centering is done using the Kernel.
%}
one_mat = ones(size(K))/size(K,1);
K_center = K - one_mat*K - K*one_mat + one_mat*K*one_mat;
clear K
%% STEP 3
%% Obtaining the low dimensional projection
% The following equation needs to be satisfied for K N*lamda*K*alpha = K*alpha 
% Thus lamda's has to be normalized by the number of points
[eigvec, eigval] = eigs(K_center);
eig_val = eigval./size(data_in,2);
%% STEP 4
% Again 1 = lamda*(alpha.alpha)
% Here '.' indicated dot product
for col = 1:size(eigvec,2)
    eigvec(:,col) = eigvec(:,col)./(sqrt(eig_val(col,col)));
end
[~, index] = sort(eig_val,'descend');
eigvec = eigvec(:,index);
eigen_val = diag(eig_val);
count = eigen_val(1,1);
new_eig_val = eigen_val(1,1);
for i = 2:size(data_in,2)
     t = (count + eigen_val(i,1));
     expl = t/sum(eigen_val);
     if expl < var/100
         new_eig_val = [new_eig_val; eigen_val(i,1)];
     
     else 
        break 
     end
     count = t;
end    
lambdas = new_eig_val;
num_dim = size(lambdas,1);
temp = zeros(num_dim,size(data_in,2));
for i = 1:num_dim
    for j = 1:size(data_in,2)
        temp(i,j) = eigvec(i,j);
    end    
end
eig_vector = real(temp);
end