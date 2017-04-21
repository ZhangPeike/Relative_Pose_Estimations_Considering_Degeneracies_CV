function Rt=Rt_from_Projection_mat_K(K,X_set,x_set)
%ZhangPeike
%Reference MVG 2nd Ed P179
%2016Äê10ÔÂ4ÈÕ15:38:07
%Input
%K: internal parameter matrix
%x_set: image of 3D points, 2D, inhomogenerous
%X_set: 3D points, inhomogenerous
%Output: Global pose parameters
% V2 2016-12-09 10:33:17
% V3 Time cost is reduced with no loop
% 2017-04-06 17:21:23
    N=size(X_set,2);
    % Normalizing
    [T2D,x_nmz]=NormalizingV2(x_set);
    [T3D,X_nmz]=Normalizing3DV2(X_set);
    X_set=[X_nmz;ones(1,N)];
    %%
    A=zeros(3*N,12);
    A(1:N,5:12)=[-X_set' repmat(x_nmz(2,:)',[1,4]).*X_set'];
    A(N+1:2*N,1:4)=X_set';
    A(N+1:2*N,9:12)=repmat((-x_nmz(1,:)'),[1,4]).*X_set';
    A(2*N+1:3*N,1:8)=[repmat((-x_nmz(2,:)'),[1,4]).*X_set' repmat(x_nmz(1,:)',[1,4]).*X_set'];
    [~,~,V]=svd(A);
    P=reshape(V(:,12),4,3)';
    P=T2D\P*T3D;
    Rt=K\P;
    ratio=nthroot(det(Rt(:,1:3)),3);
    Rt=Rt./ratio;
    %Enforce to be rotation matrix
    R=Rt(:,1:3);
    [U,~,V]=svd(R);
    Rt=[U*V',Rt(:,4)];
end