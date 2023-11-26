pic = imread('new4.jpg');
pic = double(pic)/255;
load pic_Omega
Omega = full([Ome;Ome;Ome]);
M = sparse([pic(:,:,1);pic(:,:,2);pic(:,:,3)].*Omega);

parfor i = 1
end
t = cputime;
%[X, iters, k] = SVT(M, 0.1, [0 1]);
[X, iters, k] = fastSVT_U(M, 0.1, [0 1], 50, 10);
% 将稀疏矩阵转换为完全矩阵
X_full = full(X);

% 如果矩阵的数据不是 uint8 类型，需要转换
if ~isa(X_full, 'uint8')
    % 标准化数据到 0 到 1
    X_full = X_full - min(X_full(:));
    X_full = X_full / max(X_full(:));
    % 转换到 0 到 255 范围，并转为 uint8
    X_full = uint8(X_full * 255);
end

imwrite(X_full, 'output_image.jpg');
t_SVT = cputime - t
err = sum(sum(abs([pic(:,:,1);pic(:,:,2);pic(:,:,3)]-X)))/2048/2048/3*255
