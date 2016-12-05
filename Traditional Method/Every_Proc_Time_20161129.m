%Goal:Time recording
%Read in images
clear;
K=[1520.400000 0.000000 302.320000;
    0.000000 1525.900000 246.870000;
    0.000000 0.000000 1.000000];
W=[0,-1,0;
    1,0,0;
    0,0,1];
P1=K*[1 0 0 0;...
      0 1 0 0;...
      0 0 1 0];
image{1}='G:\ComputerVision2016\3D Program\templeRing\templeR0013.png';
image{2}='G:\ComputerVision2016\3D Program\templeRing\templeR0014.png';
t0=tic;
img_1=imread(image{1});
img_2=imread(image{2});
t_elap(1)=toc(t0);
fprintf('Read in 2 images time cost ：%.1f ms\n',1000*t_elap(1));
clear image t0
t0=tic;
matches=MySIFT(img_1,img_2);
t_elap(2)=toc(t0);
fprintf('SIFT features extracting and matching in 2 images time cost ：%.1f ms\n',1000*t_elap(2));
clear t0 img_1 img_2;
t0=tic;
[F, inliers]=ransacfitfundmatrix(matches(1:2,:),matches(3:4,:),0.0015,0);
t_elap(3)=toc(t0);
fprintf('Estimating Fundamental matrix time cost ：%.1f ms\n',1000*t_elap(3));
clear t0
t0=tic;
E=K'*F*K;
t_elap(4)=toc(t0);
fprintf('Transfering Fundamental to Essential matrix：%.2f ms\n',1000*t_elap(4));
clear F t0
Rt=zeros(3,4,4);
t0=tic;
[U,D,V]=svd(E,0);
d=(D(1,1)+D(2,2))/2;
E=U*diag([d,d,0])*V';
[U,~,V]=svd(E,0);
R{1}=U*W*V';
R{2}=U*W'*V';
if det(R{1})<0
    R{1}=-R{1};
end
if det(R{2})<0
    R{2}=-R{2};    
end
t{1}=U(:,3);
t{2}=-t{1};
l=1;
for j=1:2
    for k=1:2
        Rt(:,:,l)=[R{j} t{k}];
        l=l+1;
    end
end
t_elap(5)=toc(t0);
fprintf('Getting 4 poses from Essential matrix：%.1f ms\n',1000*t_elap(5));
clear d t0 D
t0=tic;
for m=4:-1:1
    P2=K*Rt(:,:,m);
    Vote(m)=0;
    %注意函数的输入，三维点坐标    
    X=Triangulation(P1,P2,matches(:,inliers));
    %化为非齐次
    X = X(1:3,:) ./ X([4 4 4],:);
    for n=1:size(X,2)
        mmcos(n)=Rt(3,1:3,m)*(X(:,n)+Rt(:,1:3,m)'*Rt(:,4,m));        
        Vote(m)=(X(3,n)>0 & mmcos(n)>0)+Vote(m);
    end    
end
[~,index]=max(Vote);
RightRt=Rt(:,:,index); 
t_elap(6)=toc(t0);
fprintf('Getting the right pose from 4pose  ：%.1f ms\n',1000*t_elap(6));
for i=1:4
    RtN{i}=Rt(:,:,i);
end
clear t0 index i j k
t0=tic;
npts=size(matches(:,inliers),2);
x1=[matches(1:2,inliers);ones(1,npts)];
x2=[matches(3:4,inliers);ones(1,npts)];
x1n=K\x1;
x2n=K\x2;
Rt_NM = RT_Check( x1n,x2n,RtN );
t_elap(7)=toc(t0);
fprintf('New method getting the right pose from 4 pose  ：%.1f ms\n',1000*t_elap(7));