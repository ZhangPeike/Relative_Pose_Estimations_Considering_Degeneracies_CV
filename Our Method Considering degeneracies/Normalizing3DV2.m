function [ T,Tpoint ] = Normalizing3DV2( point )
%对3D点云的归一化处理
%Normalization
%point是非齐次坐标
%过程，首先使图像点的中心平移到坐标原点，
%平移和缩放用3*3矩阵表示
%参考MVG P109
%先求中心，再求平均距离
%T=[*,*,*;*,*,*;0,0,1]
%作者：张培科
%2016年 2月6日10:32:25 
%2016年10月4日10:59:44
%Version 3
%Inhomogeneous points
%V2 2017-04-06 15:50
%%
%求点集的中心坐标
centroid=mean(point);
newpoint(1,:)=point(1,:)-centroid(1);
newpoint(2,:)=point(2,:)-centroid(2);
newpoint(3,:)=point(3,:)-centroid(3);
dist=sqrt(newpoint(1,:).^2+newpoint(2,:).^2+newpoint(3,:).^2);
meandist=mean(dist(:));
%缩放点集坐标，以使平均长度为根号3         
scale=sqrt(3)/meandist;
T=[scale,0,0,-scale*centroid(1);
   0,scale,0,-scale*centroid(2);
   0,0,scale,-scale*centroid(3);
   0,0,0,1];
%%
Newpoint=T*[point;ones(1,size(point,2))];
%归一化后的点的非齐次坐标
Tpoint=Newpoint(1:3,:);
end 