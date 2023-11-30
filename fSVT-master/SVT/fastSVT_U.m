function [X, iters, k] = fastSVT_U(M, tol, ran, i_reuse, q_reuse, delta)
% The ran is a new parameter for control the result matrix value's range
iters = 0
k = 0
[m, n]= size(M);
Omega = spones(M);
Ns = sum(sum(Omega))
if nargin == 5
    delta = 1.2*m*n/Ns;
end
[xi, yi, ~] = find(Omega);
tau=5*n;
l=5;
i_max=300;
PM = M;
normPM2= svds(PM, 1);
normPM= norm(PM, 'fro');
k0=ceil(tau/delta/normPM2);
Y0 = k0*delta*PM;

dec = 0;
r = 0;
p = 2;
q = 0;
err_before = 1000;

for i = 1:i_max
    r_before = r;
    r=r+1;
    if mod(i,50) == 0
        delta = delta/1.1;    % This is very important for stablizing convergence.
    end                       % Don't know the reas
    if i > i_reuse && q < q_reuse
        [U,S,V,Q] = rSVDBKIr(Y0, r, p, U);
        q = q + 1;
    else
        [U,S,V,Q] = rSVDBKIr(Y0, r, p);
        q = 0;
    end
    while S(1)>tau
        r=r+l;
        [U,S,V,Q] = rSVDBKIr(Y0, r, p);
    end
    for j= 1:r,
        if S(j)> tau,
            break
        end
    end
    r_max = r;
    r = max(r_max-j+1, r_before);
    x = r_max-r+1:r_max;
    S(x) = S(x)-ones(r,1)*tau;
    parfor j = x
       U(:, j) = U(:, j)*S(j);
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
    X = sparse(xi, yi, x_now, m, n);
    PX= X-PM;
    err = norm(PX, 'fro')/normPM;
    if err <= tol,
        X = U(:, x)*V(:, x)';
        X(X<ran(1)) = ran(1);
        X(X>ran(2)) = ran(2);
        k = r;
        iters = i;
        break;
    end
        if err > err_before
        dec = 0;
        p = p + 1;
        q = 10;
    else
        if p <= 5
            dec = 0;
        else
            dec = dec + 1;
            if dec == 10
                p = p - 1;
                dec = 0;
            end
        end
        end
    err_before = err;
    disp([i, r, err, p]);
    Y0= Y0 - delta*PX;
    t = cputime;
%[X, iters, k] = SVT(M, 0.1, [0 1]);
%[X, iters, k] = fastSVT_U(M, 0.1, [0 1], 50, 10);
% 将稀疏矩阵转换为完全矩阵
if mod(i,50)==0
X_full = full(X);

len = 2048
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
filename = sprintf('output_image_fast_SVT_U_i%d.jpg', i);
imwrite(image, filename);
t_SVT = cputime - t
disp(t_SVT)
end
end
end