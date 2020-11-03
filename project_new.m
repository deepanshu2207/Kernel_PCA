%% STEP 5
function new_projected_data = project_new(x_new, data_in, type, par1, par2, alphas)
%{
x_new = New data to be projected on the transformed data
data_in = Data used for previous function
type = Choose the kernel function to use: 'G' for Gaussian, 'P' for Polynomial, 'S' for Sigmoid
par1, par2 = values of parameters for the chosen kernel
alphas = Output of previous function (k selected eigenvectors)
new_projected_data (OUTPUT) =    
%}    
   % Normalization of new data
   old_variance = var(data_in');
   old_mean = mean(data_in');
  for i=1:size(x_new,2)
   x_new(:,i) = (x_new(:,i)-old_mean')./sqrt(old_variance');
  end 
   % Final Step
  temp = zeros(size(data_in,2),size(x_new,2));
  temp2 = zeros(size(alphas,1),size(alphas,2));
  new_projected_data = zeros(size(alphas,1),size(x_new,2));
    if type == 'G'
      for row = 1:size(data_in,2)
          for col = 1:size(x_new,2)
              temp(row,col) = sum(((data_in(:,row) - x_new(:,col)).^2));
          end    
      end
      for i = 1: size(data_in,2)
          for j = 1:size(x_new,2)
       temp2(: , i) = (temp(i,j).*alphas(:,i));
       new_projected_data(:,j) = sum(temp2,2);
          end
      end    
   end
   
   if type == 'P'
      for row = 1:size(data_in,2)
           for col = 1:size(x_new,2)
               temp(row,col) = (dot(data_in(:,row), x_new(:,col))+ par2)^par1;
           end
      end
       for i = 1: size(data_in,2)
             for j = 1:size(x_new,2)
                 temp2(: , i) = (temp(i,j).*alphas(:,i));
                 new_projected_data(:,j) = sum(temp2,2);
             end    
      end 
   end
   
    if type == 'S'
       for row = 1:size(data_in,2)
            for col = 1:size(x_new,2)
                temp(row,col) = tanh(par1*(data_in(:,row)'*x_new(:,col)) + par2);
            end    
       end
         for i = 1: size(data_in,2)
             for j = 1:size(x_new,2)
                 temp2(: , i) = (temp(i,j).*alphas(:,i));
                 new_projected_data(:,j) = sum(temp2,2);
             end    
        end 
    end
end