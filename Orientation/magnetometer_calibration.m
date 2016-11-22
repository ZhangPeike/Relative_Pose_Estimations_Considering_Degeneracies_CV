function [D d] = magnetometer_calibration(y)
% S1: Initial estimate for both estimations
% S2: Algorithm 1 (Norm matching)
% S3: Algorithm 2 (MLE)

%% S1: Initial estimate for both estimations
n = length(y);
Y = kron(y(:,1)',y(:,1)');
for i = 2:n
    Y = [Y; kron(y(:,i)',y(:,i)')];
end
Y(:,[4 7 8]) = Y(:,[2 3 6]) + Y(:,[4 7 8]); % As A is symmetric, eliminate lower-left 3 redundant elements and estimate upper-right elements only.
Y(:,[2 3 6]) = [];
Y = [Y y' ones(n,1)];
[u s v] = svd(Y); % SVD: Y = u*s*v'
tmeta = v(:,end); % vn: last eigenvector of V
A = [[tmeta(1) tmeta(2) tmeta(4)]' [tmeta(2:3); tmeta(5)] tmeta(4:6)];
b = tmeta(7:9); c = tmeta(10);

alfa = 1/(0.25*b'*inv(A)*b - c);
A = alfa*A; b = alfa*b; c = alfa*c;
h0 = -0.5*inv(A)*b; R0 = chol(A); % A = R'*R
% % %% norm checking
% x = R0*(y-repmat(h0,1,n));
% figure,
% [s1 s2 s3] = sphere; surf(s1,s2,s3); grid
% hold on; plot3(x(1,:),x(2,:),x(3,:),'b.'); hold off; grid;
% alpha(0.5); % make transparent
% axis equal;
%
% figure, plot(vmmag(x),'b'); title('Points Norm (by initial estimate)')

alfa = 0.;
R0 = R0.*(1+alfa*sign(randn(3)));
h0 = h0.*(1+alfa*sign(randn(3,1)));

%% S2: Algorithm 1: Norm matching
disp('## Algorithm 1: Norm matching');
% %% Iterative estimate by Newton method
T = R0;
x = reshape(T,9,1); x([2 3 6]) = [];
x = [x; h0]; % x0

clear obj;
disp('obj(x) = ');
for k = 1:5
    tmp = func(x);
    obj(k) = sum(sum(tmp.*tmp,1)); % objective function value at current x
    disp(['          ' num2str(obj(k)) ': ' num2str(x')]);

    T = [[x(1) 0 0]' [x(2:3); 0] x(4:6)]; d = x(7:9);
    J = 0; H = 0;
    for i = 1:n
        u = y(:,i) - d;
        Jt = 4*(norm(T*u)^2-1)*kron(u,T*u);
        Jd = -4*(norm(T*u)^2-1)*T'*T*u;
        J = J + [Jt; Jd];
        Htt = 8*kron(u*u',T*u*u'*T') + 4*(norm(T*u)^2-1)*kron(u*u',eye(3));
        Htd = -8*kron(u,T*u)*u'*T'*T - 4*(norm(T*u)^2-1)*(kron(eye(3),T*u) + kron(u,T));
        Hdd = 4*(norm(T*u)^2-1)*T'*T + 8*T'*T*u*u'*T'*T;
        H = H + [Htt Htd; Htd' Hdd];
    end
    % remove elements for 3 left-lower elements
    J([2 3 6]) = [];
    H(:,[2 3 6]) = []; H([2 3 6],:) = [];

    dx = -H\J;%-inv(H)*J;
    x = x + dx;
end
T = [[x(1) 0 0]' [x(2:3); 0] x(4:6)]; d = x(7:9);
tmp = func(x);
obj(end+1) = sum(sum(tmp.*tmp,1)); % objective function value at current x
disp(['          ' num2str(obj(k)) ': ' num2str(x')]);

% % Monte Carlo recording
% fprintf(fh, '%f   ', [x; obj']);  fprintf(fh,'\n');

% %% norm checking
x = T*(y-repmat(d,1,length(y)));
mag = vmmag(x);
figure,
subplot(2,1,1); plot(mag); title('Points Norm (by NM)')
subplot(2,1,2); histfit(mag-1);
disp(['mean: ' num2str(mean(mag-1)) '  std: ' num2str(std(mag))]);
disp(['Calibrated parameters:']);
T
d
D = T;

% %% S3: Algorithm 2: MLE
% disp('## Algorithm 2: MLE');
% % %% norm checking
% m = R0*(y-repmat(h0,1,n));
% % %% Iterative estimate by Newton method
% D = inv(R0);
% x = reshape(D,9,1); x([2 3 6]) = [];
% lamda = zeros(n,1);
% x = [x; h0; reshape(m,3*n,1); lamda]; % dim = 4*n+9
% x0 = x;
% 
% clear obj;
% disp('obj(x) = ');
% for k = 1:5
%     D = [[x(1) 0 0]' [x(2:3); 0] x(4:6)]; d = x(7:9);
%     m = reshape(x(9+1:9+3*n),3,n); lamda = x(9+3*n+1:end);
%     tmp = y-D*m-repmat(d,1,n);
%     obj(k) = sum(sum(tmp.*tmp,1)); % objective function value at current x
%     disp(['          ' num2str(obj(k)) ': ' num2str(x(1:9)')]);
% 
%     J = zeros(4*n+9,1); H = zeros(4*n+9);
%     for i = 1:n
%         u = y(:,i) - d;
%         J(1:9) = J(1:9) - 2*kron(m(:,i),u-D*m(:,i)); % JD
%         J(10:12) = J(10:12) - 2*(u-D*m(:,i)); % Jd
%         J(12 + ((i-1)*3+1:i*3)) = -2*D'*(u-D*m(:,i)) + 2*lamda(i)*m(:,i); % Jm_k
%         J(12+3*n + i) = m(:,i)'*m(:,i) - 1; % Jlamda_k
% 
%         H(1:9,1:9) = H(1:9,1:9) + 2*kron(m(:,i)*m(:,i)',eye(3)); % JDD
%         H(1:9,10:12) = H(1:9,10:12) + 2*kron(m(:,i),eye(3)); % JDd
%         H(1:9,12 + ((i-1)*3+1:i*3)) = 2*(kron(m(:,i),eye(3))*D - kron(eye(3),u-D*m(:,i))); % JDm_k
% 
%         H(10:12,1:9) = H(1:9,10:12)'; % JdD
%         H(10:12,10:12) = H(10:12,10:12) + 2*eye(3); % Jdd
%         H(10:12,12 + ((i-1)*3+1:i*3)) = 2*D; % Jdm_k
% 
%         H(12 + ((i-1)*3+1:i*3),1:9) = H(1:9,12 + ((i-1)*3+1:i*3))'; % Jm_k/D
%         H(12 + ((i-1)*3+1:i*3),10:12) = H(10:12,12 + ((i-1)*3+1:i*3))'; % Jm_k/d
%         H(12 + ((i-1)*3+1:i*3),12 + ((i-1)*3+1:i*3)) = 2*(D'*D + lamda(i)*eye(3)); % Jm_k^2
%         H(12 + ((i-1)*3+1:i*3),12+3*n + i) = 2*m(:,i); % Jm_k/lamda_k
% 
%         H(12+3*n + i,12 + ((i-1)*3+1:i*3)) = H(12 + ((i-1)*3+1:i*3),12+3*n + i)'; %Jlamda_k/m_k
%     end
%     % remove elements for 3 left-lower elements
%     J([2 3 6]) = [];
%     H(:,[2 3 6]) = []; H([2 3 6],:) = [];
% 
%     dx = -H\J;%-inv(H)*J;
%     x = x + dx;
% end
% D = [[x(1) 0 0]' [x(2:3); 0] x(4:6)]; d = x(7:9);
% m = reshape(x(9+1:9+3*n),3,n); lamda = x(9+3*n+1:end);
% tmp = y-D*m-repmat(d,1,n);
% obj(end+1) = sum(sum(tmp.*tmp,1)); % objective function value at current x
% disp(['          ' num2str(obj(k)) ': ' num2str(x(1:9)')]);
% 
% D = inv(D);
% tmp = reshape(D,9,1); tmp([2 3 6]) = [];
% x(1:6) = tmp;
% % % Monte Carlo recording
% % fprintf(fh, '%f   ', [x(1:9); obj']);  fprintf(fh,'\n');
% 
% % %% norm checking
% x = D*(y-repmat(d,1,length(y)));
% mag = vmmag(x);
% figure,
% subplot(2,1,1); plot(mag); title('Points Norm (by ML)')
% subplot(2,1,2); histfit(mag-1);
% disp(['mean: ' num2str(mean(mag-1)) '  std: ' num2str(std(mag))]);
% 
% disp(['Calibrated parameters:']);
% D
% d


