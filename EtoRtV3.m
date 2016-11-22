function RightRt = EtoRtV3(E,K1,K2,matches)
%Recover the relative pose from the essential matrix
%Input：E，essential matrix,3*3 
%K1，camera 1 internal parameter matrix
%K2，camera 2 internal parameter matrix
%matches，4*n for n point correspondence
%Output: the relative pose [R t]
%Reference Multiple View Geometry in Computer Vision, 2nd edition,P257
%Author: Peike Zhang, Yuanxin Wu, Qi Cai and Danping Zou
%2016
%Version 2
%Version：3
W=[0,-1,0;
    1,0,0;
    0,0,1];
Z=[0,1,0;
  -1,0,0;
   0,0,0];
[U,w,V]=svd(E);
w_diag=(w(1,1)+w(2,2))/2;
%2016-08-12 11:03:48
w(1,1)=w_diag;
w(2,2)=w_diag;
w(3,3)=0;
U=det(U)*U;
V=det(V)*V;
E=U*w*V';
[U,~,V]=svd(E);
%MVG P258 Equation9.14
S=U*Z*U';
%U=det(U)*U;
%V=det(V)*V;
R{1}=U*W*V';
R{2}=U*W'*V';
%Nister 5points algo Theorem 2
%{
disp('E2Rt Theorem2 by Nister');
fprintf('det(U)%.3f\n',det(U));
fprintf('det(V)%.3f\n',det(V));
fprintf('det(R1)%.3f\n',det(R{1}));
fprintf('det(R2)%.3f\n',det(R{2}));
%}
%增加的行列式判断，增加后结果正确
if det(R{1})<0
    R{1}=-R{1};
end
if det(R{2})<0
    R{2}=-R{2};
end
t{1}=U(:,3);
t{2}=-t{1};
%Find out the right pose by the positive depth constraint.
l=1;
for j=1:2
    for k=1:2
        Rt(:,:,l)=[R{j} t{k}];
        l=l+1;
    end
end
%准备试算
P1=K1*[1 0 0 0;...
     0 1 0 0;...
     0 0 1 0];
for i=1:4
    P2=K2*Rt(:,:,i);
    Vote(i)=0;
    %注意函数的输入，三维点坐标    
    X=Triangulation(P1,P2,matches);
    %化为非齐次
    X = X(1:3,:) ./ X([4 4 4],:);
    for m=1:size(X,2)
        %齐次坐标转化为非齐次
        %X(:,m)=X(1:3,m)./X([4 4 4],m);
        %判断从第二个相机焦点到X点的矢量(X(:,m)-Rt(:,4,i))与第二个相机的光轴的点积
        %R矩阵的第三行Rt(3,1:3,i)即Z轴方向的单位矢量在世界坐标系的坐标
        %四个解的判断方法与SFMedu的比较理解
        %2016年5月25日17:00:04，更正，第二个相机的角点坐标C=-R'*t
        %可增加 重构点与两个光轴的夹角
        mmcos(m)=Rt(3,1:3,i)*(X(:,m)+Rt(:,1:3,i)'*Rt(:,4,i));        
        Vote(i)=(X(3,m)>0&mmcos(m)>0)+Vote(i);
    end    
end
%disp('备选解');
%disp(Vote);
%disp(dcm2eul(R{1})*180/pi);
%disp(dcm2eul(R{2})*180/pi);
[~,index]=max(Vote);
RightRt=Rt(:,:,index);
end