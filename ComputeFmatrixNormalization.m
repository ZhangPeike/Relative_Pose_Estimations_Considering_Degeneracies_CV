function [ F ] = ComputeFmatrixNormalization(matches)
%Computing Fundamental Matrix
%Peike Zhang, Yuxin Wu, Qi Cai, Danping Zou
%Reference Multiple View Geometry in Computer Vision Chapter 13
[a,b]=size(matches);
x1=matches(1:2,:);
x2=matches(3:4,:);
%归一化，归一化函数有问题，函数结果无转换后的点
[T1,x1]=MyNormalizing(x1);
[T2,x2]=MyNormalizing(x2);
%
%{
 A = [x2(1,:)'.*x1(1,:)'   x2(1,:)'.*x1(2,:)'  x2(1,:)' ...
         x2(2,:)'.*x1(1,:)'   x2(2,:)'.*x1(2,:)'  x2(2,:)' ...
         x1(1,:)'             x1(2,:)'            ones(npts,1) ];       

%}

A=[x2(1,1)*x1(1,1),x2(1,1)*x1(2,1),x2(1,1),x2(2,1)*x1(1,1),x2(2,1)*x1(2,1),x2(2,1),x1(1,1),x1(2,1),1];
for i=1:b-1
    A=[A;x2(1,1+i)*x1(1,1+i),x2(1,1+i)*x1(2,1+i),x2(1,1+i),x2(2,1+i)*x1(1,1+i),x2(2,1+i)*x1(2,1+i),x2(2,1+i),x1(1,1+i),x1(2,1+i),1];
end
%SVD
[U,S,V]=svd(A);
%F = reshape(V(:,9),3,3)';
tempF=[V(1,9),V(2,9),V(3,9);
   V(4,9),V(5,9),V(6,9);
   V(7,9),V(8,9),V(9,9)];
%rank(F)=2
[U1,S1,V1]=svd(tempF);
tempF=U1*diag([S1(1,1),S1(2,2),0])*V1';
F=T2'*tempF*T1;
