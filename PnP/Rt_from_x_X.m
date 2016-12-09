function Rt=Rt_from_x_X(K,X_set,x_set)
%ZhangPeike
%2016年10月4日15:38:07
%Input
%K: internal parameter matrix
%x_set: image of 3D points, 2D, inhomogenerous
%X_set: 3D points, inhomogenerous
%Output: Global pose parameters
% Normalizing
N=size(X_set,2);
[T2D,x_nmz]=MyNormalizing(x_set);
[T3D,X_nmz]=MyNormalizing3D(X_set);
X_set=[X_nmz;ones(1,N)];
A=[zeros(1,4),-X_set(:,1)',x_nmz(2,1)*X_set(:,1)';
   X_set(:,1)',zeros(1,4),-x_nmz(1,1)*X_set(:,1)';
    -x_nmz(2,1)*X_set(:,1)',x_nmz(1,1)*X_set(:,1)',zeros(1,4)];
for i=2:N
    A=[A;
      zeros(1,4),-X_set(:,i)',x_nmz(2,i)*X_set(:,i)';
      X_set(:,i)',zeros(1,4),-x_nmz(1,i)*X_set(:,i)';
      -x_nmz(2,i)*X_set(:,i)',x_nmz(1,i)*X_set(:,i)',zeros(1,4)
      ];   
    %2016年9月27日16:09:06
    %2016年10月4日11:32:20
    [~,~,V]=svd(A);
    P=[V(1,12),V(2,12),V(3,12),V(4,12);
       V(5,12),V(6,12),V(7,12),V(8,12);
       V(9,12),V(10,12),V(11,12),V(12,12);];
    P=T2D\P*T3D;
    Rt=K\P;
end
end