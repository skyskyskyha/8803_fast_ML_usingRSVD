function [X, iters, k] = SVT(M, tol, ran, delta)
% This is the method in [Cai et al., 2010] using svds
% The ran is a new parameter for control the result matrix value's range
[m, n]= size(M);
Omega = spones(M);
Ns = sum(sum(Omega))
if nargin == 3
    delta = 1.2*m*n/Ns;
end
[xi, yi, ~] = find(Omega);
tau=5*n;
l=5;
i_max=1000;
PM = M;
t=cputime;
normPM2= svds(PM, 1);
normPM= norm(PM, 'fro');
k0=ceil(tau/delta/normPM2);
Y0 = k0*delta*PM;

r=0;
for i = 1:i_max
    r_before = r;
    r=r+1;
    if mod(i,50) == 0
        delta = delta/1.1;    % This is very important for stablizing convergence.
    end                       % Don't know the reas
    [U, S, V] = svds(Y0, r);
    while S(r,r)>tau
        r=r+l;
        [U,S,V]=svds(Y0, r);
    end
    for j= r:-1:1,
        if S(j, j)> tau,
            break
        end
    end
    r = max(j, r_before);
    x = 1:r;
    s = diag(S(x, x)) - ones(r, 1)*tau;
    parfor j = x
        U(:, j) = U(:, j)*s(j);
    end
    x_now = zeros(Ns, 1);
    parfor j = 1:Ns
        temp = U(xi(j), x)*V(yi(j), x)';
        if temp < ran(1)
            x_now(j) = ran(1);
        elseif temp > ran(2)
            x_now(j) = ran(2);
        else
            x_now(j) = temp;
        end
    end
    X = sparse(xi, yi, x_now);
    PX= X-PM;
    err= norm(PX, 'fro')/normPM;
    if err <= tol,
        X = U(:, x)*V(:, x)';
        X(X<ran(1)) = ran(1);
        X(X>ran(2)) = ran(2);
        k = r;
        iters = i;
        break;
    end
    disp([i, r, err]);
    Y0= Y0 - delta*PX;
if mod(i,50)==0
X_full = full(X);

len = 100
% 分离三个颜色通道
R = X_full(1:len*1, :);         % 红色通道
G = X_full(len+1:len*2, :);      % 绿色通道
B = X_full(len*2+1:end, :);       % 蓝色通道


% 重组成一个三维彩色图像
image = cat(3, R, G, B);

% 确保图像数据类型是 uint8
if ~isa(image, 'uint8')
    % 标准化到 0 到 1 范围
    image = image - min(image(:));
    image = image / max(image(:));
    % 转换到 0 到 255 范围，并转为 uint8
    image = uint8(image * 255);
end

% 使用 imwrite 保存图像
filename = sprintf('output_image_common_SVT_i%d.jpg', i);
imwrite(image, filename);
t_SVT = cputime - t
disp(t_SVT)
end
end
end