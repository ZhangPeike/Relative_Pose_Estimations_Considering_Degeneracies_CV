function  XSet  = Triangulation(P1,P2,matches)
%基本思想，构造关于X坐标的方程
%最优化问题
%注意，函数输入matches，4*n，前两行对应P1，后两行对应P2。
%输出：齐次坐标
%Author:Peike Zhang
%x叉乘PX=0，构造方程，MVG P312
for i=size(matches,2):-1:1
A=[matches(1,i)*P1(3,:)-P1(1,:);
   matches(2,i)*P1(3,:)-P1(2,:);
   matches(1,i)*P1(2,:)-matches(2,i)*P1(1,:);
   matches(3,i)*P2(3,:)-P2(1,:);
   matches(4,i)*P2(3,:)-P2(2,:);
   matches(3,i)*P2(2,:)-matches(4,i)*P2(1,:)];
[~,~,V]=svd(A,0);
%V的最后一列为解集，每4个为三维坐标
XSet(:,i)=V(:,end);
end
end

