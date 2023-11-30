len = 100
pic = imread('img7.png');
pic = double(pic)/255;
load pic_Omega
Omega = full([Ome;Ome;Ome]);
M = sparse([pic(:,:,1);pic(:,:,2);pic(:,:,3)].*Omega(1:len*3, 1:len));

parfor i = 1
end
t = cputime;
[X, iters, k] = SVT(M, 0.1, [0 1]);
%[X, iters, k] = fastSVT_Q(M, 0.2, [0 1], 50, 10);
% 将稀疏矩阵转换为完全矩阵
X_full = full(X);


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
imwrite(image, 'output_image_small_SVT.jpg');
t_SVT = cputime - t
err = sum(sum(abs([pic(:,:,1);pic(:,:,2);pic(:,:,3)]-X)))/2048/2048/3*255
