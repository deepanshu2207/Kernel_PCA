# Kernel_PCA
This was my college project course in which I learnt about Kernel PCA and its application in Chemical Engineering.

There are two matlab files, namely, kernelpca.m and project_new.m. Both of the codes have been properly commented as per the algorithm written in project report.

1. User need to choose the type and parameters of the kernel function. Now, there are 3 choices of the kernel function in code, that is, Gaussian, Sigmoid, and Polynomial.
2. Instead of mentioning the number of PCs, user need to specify the amount of variance (in percentage) he/she wants to retain in the transformed data.
3. Data is standardized before doing any calculations.
4. I have made 2 functions, one to apply kernel PCA on given data (kernelpca.m) and other (project_new.m) to project any new data on the transformed data. So, if we have new data to project, then we need to use both the functions, otherwise, we just apply the first function for kernel PCA.
5. As we are giving the input as variance we want to retain, we now know how many eigenvectors are required for the given amount of variance which is to be retained.
