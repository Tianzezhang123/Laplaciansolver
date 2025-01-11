e_2 = zeros(10,1);
times = zeros(10,1);
for k = 1:10
    N = 10*k;
    
    % 在[-1,1]x[0,1]上构建网格
    x = linspace(-1, 1, 2*N+1);
    y = linspace(-1, 1, 2*N+1);
    % 定义底边函数
    BottomBoundary = @(x) 1/2*(cos(x*pi)-1);
    % 每列的格点个数
    ynumb = zeros(1,2*N+1,'uint64');
    for i = 1:2*N+1
        if i ~= N+1
            ynumb(i) = 1+floor(N*(1-max(-1,BottomBoundary(x(i)))));
        else
            ynumb(i) = N+1;
        end
    end
    totnum = sum(ynumb);
    % 边值函数（先取一阶近似）
    boundaryFunc = @(x,y) x^2-y^2;
    
    % 真解
    ureal = zeros(2*N+1,2*N+1);
    for i = 1:2*N+1
        for j = 1:2*N+1-ynumb(i)
            ureal(i,j) = nan;
        end
        for j = 2*N+2-ynumb(i):2*N+1
            ureal(i,j) = boundaryFunc(x(i),y(j));
        end
    end
    
    % 数值解
    tic()
    u = sin2dLap(N, boundaryFunc);
    times(k) = toc();
    
    % 可视化解
    % [X, Y] = meshgrid(x, y);
    % surf(X, Y, u');
    % title('Real Solution to Laplace Equation in Domain with Spine');
    % xlabel('x');
    % ylabel('y');
    % zlabel('u');
    
    % 计算误差
    e_2k = 0;
    ureal_2=0;
    for i = 1:2*N+1
        for j = 1:2*N+1
            if isnan(u(i,j)) || isnan(ureal(i,j))
            else
                e_2k = e_2k + (u(i,j)-ureal(i,j))^2;
                ureal_2 = ureal_2+ureal(i,j)^2;
            end
        end
    end
    e_2k = sqrt(e_2k/ureal_2);
    e_2(k) = e_2k; 
end

disp(times)

x = 20*(1:10); 
y = e_2;  

loglog(x, y, 'm*-');
grid on;
title('Error of the 4-points method with sin boundary');
xlabel('mesh number (log scale)');
ylabel('eror (log scale)');