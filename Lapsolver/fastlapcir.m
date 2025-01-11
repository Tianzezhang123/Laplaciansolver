e_2 = zeros(10,1);
times = zeros(10,1);
for n = 1:10
    N = 100*n; % 边界点数量
    theta = linspace(0, 2*pi, N+1); theta(end) = []; % 边界参数化
    x_boundary = cos(theta); % 边界点的 x 坐标
    y_boundary = sin(theta); % 边界点的 y 坐标
    
    % Dirichlet 边界条件
    f = @(x, y) x.^2 - y.^2; % 边界条件
    f_boundary = f(x_boundary, y_boundary); % 边界上的值
    
    % 法向量 (单位圆上的外法向)
    nx = x_boundary; 
    ny = y_boundary;
    
    % 双层势核函数的法向导数
    kernel = @(x, y, nx, ny, tx, ty) ...
        ((nx .* (x - tx) + ny .* (y - ty)) ./ ((x - tx).^2 + (y - ty).^2)) / (2*pi);
    
    tic()
    % 构造积分方程矩阵
    A = zeros(N, N);
    for i = 1:N
        for j = 1:N
            if i ~= j
                A(i, j) = kernel(x_boundary(i), y_boundary(i), ...
                                 nx(i), ny(i), ...
                                 x_boundary(j), y_boundary(j))*2*pi/N;
            else
                A(i, j) = 1/2*(1+1/N); % 自身积分项
            end
        end
    end
    
    % 求解密度函数 sigma
    sigma = A \ f_boundary(:);
    
    % 计算单位圆内部的解
    n_r = 10; % 半径分布
    n_theta = 50; % 角度分布
    r = linspace(0, 1-1/n_r, n_r);
    theta_grid = linspace(0, 2*pi, n_theta);
    [R, Theta] = meshgrid(r, theta_grid);
    X = R .* cos(Theta);
    Y = R .* sin(Theta);
    
    % 计算解
    U = zeros(size(X));
    for i = 1:numel(X)
        for j = 1:N
            U(i) = U(i) - kernel(X(i), Y(i), x_boundary(j), y_boundary(j), ...
                                 x_boundary(j), y_boundary(j)) * sigma(j) *2*pi/N;
        end
    end
    times(n)= toc();
    
    % 真解
    U_true = f(X, Y);
    
    % 2-norm of error
    e_2(n) = norm(U - U_true)/norm(U_true);
end

x = 100*(1:10); 
y = e_2;  

loglog(x, y, 'm*-'); 
grid on;
title('Error of the integral approach with circular boundary');
xlabel('mesh number (log scale)');
ylabel('2-norm eror (log scale)');
