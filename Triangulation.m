function  XSet  = Triangulation(P1,P2,matches)
%Input:projection matrix,P1,P2,both 3*4;matches,4*n for n point
%correspondence
%Author: Peike Zhang, Yuanxin Wu, Qi Cai and Danping Zou
%2016
%Reference Multiple View Geometry in Computer Vision, 2nd edition,P312
for i=1:size(matches,2)
A=[matches(1,i)*P1(3,:)-P1(1,:);
   matches(2,i)*P1(3,:)-P1(2,:);
   matches(1,i)*P1(2,:)-matches(2,i)*P1(1,:);
   matches(3,i)*P2(3,:)-P2(1,:);
   matches(4,i)*P2(3,:)-P2(2,:);
   matches(3,i)*P2(2,:)-matches(4,i)*P2(1,:)];
[~,~,V]=svd(A,0);
XSet(:,i)=V(:,end);
end
end

