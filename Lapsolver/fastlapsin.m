e_2 = zeros(10,1);
times = zeros(10,1);
for n = 1:10
    halfN = 100*n; % mesh number is 4*N
    N=2*halfN;
    
    % 定义底边函数
    BottomBoundary = @(x) 1/2*(cos(x*pi)-1);
    
    % first and second derivitive
    bounddif1 = @(x) -pi/2*(sin(x*pi));
    bounddif2 = @(x) -pi^2/2*(cos(x*pi));
    
    % compute the curvature
    cur =@(x) bounddif2(x)/ (1 + (bounddif1(x))^2)^(3/2);
    
    % the boundary points
    pts = zeros(4*N,2); 
    ny = zeros(4*N,2);
    curs = zeros(4*N,1);
    ds = ones(4*N,1)*2/N;
    
    pts(1:N,1) = ones(N,1);
    pts(N+1:2*N,1) = linspace(1-2/N,-1,N);
    pts(2*N+1:3*N,1) = -1*ones(N,1);
    pts(3*N+1:4*N,1) = linspace(-1+2/N,1,N);
    
    pts(1:N,2) = linspace(-1+2/N,1,N);
    pts(N+1:2*N,2) = ones(N,1);
    pts(2*N+1:3*N,2) = linspace(1-2/N,-1,N);
    for i = 3*N+1:4*N
        botvalue = BottomBoundary(pts(i,1));
        pts(i,2) = botvalue;
        ny(i,:) = [bounddif1(pts(i,1)),-1];
        normny = 1+bounddif1(pts(i,1))^2;
        ny(i,:) = ny(i,:)/sqrt(normny);
        curs(i) = cur(pts(i,1))/N/8;
        ds(i) = ds(i)*sqrt(normny);
    end
    
    for s = 1:4
        curs(N*s) = 1/4;
    end
    
    % the normal derivitive
    ny(1:N,1) = ones(N,1);
    ny(N+1:2*N,1) = zeros(N,1);
    ny(2*N+1:3*N,1) = -1*ones(N,1);
    
    ny(1:N,2) = zeros(N,1);
    ny(N+1:2*N,2) = ones(N,1);
    ny(2*N+1:3*N,2) = zeros(N,1);
    
    % boundary function
    boundaryFunc = @(x,y) x^2-y^2;
    
    % coefficients in the integral equation
    tic()
    alpha = (1/2)*ones(size(pts,1),1)+curs;
    dygreen = zeros(size(pts,1),size(pts,1));
    for i = 1:size(pts,1)
        for j = 1:size(pts,1)
            if i ~= j
                dygreen(i,j) = -lapdygreen(pts(i,:),pts(j,:),ny(j,:))*ds(j);
            else
                dygreen(i,j) = 0;
            end
        end
    end
    kernel = dygreen;
    
    f = zeros(size(pts,1),1);
    for i = 1:size(pts,1)
        f(i) = boundaryFunc(pts(i,1),pts(i,2));
    end
    
    mu = solveieqn(kernel,alpha,f);
    
    % mesh grid of solution
    M=10;
    x = linspace(-1,1,2*M+1);
    y = linspace(-1,1,2*M+1);
    u = zeros(2*M+1,2*M+1);
    
    % grid point in each colum
    ynumb = zeros(1,2*M+1,'uint64');
    for i = 1:2*M+1
        if i ~= M+1
            ynumb(i) = 1+floor(M*(1-max(-1,BottomBoundary(x(i)))));
        else
            ynumb(i) = M+1;
        end
    end
    totnum = sum(ynumb);
    
    % exact solution
    ureal = zeros(2*M+1,2*M+1);
    for i = 1:2*M+1
        for j = 1:2*M+1-ynumb(i)
            ureal(i,j) = nan;
        end
        for j = 2*M+2-ynumb(i):2*M+1
            ureal(i,j) = boundaryFunc(x(i),y(j));
        end
    end
    ureal = real(ureal);
    
    % Integrate to get the solution
    muldygreen = zeros(size(pts,1),1);
    for i = 2:2*M
        for j = 1:2*M+1-ynumb(i)
            ureal(i,j) = nan;
        end
        for j = 2*M+2-ynumb(i):2*M
            point = [x(i),y(j)];
            for k = 1:size(pts,1)
                muldygreen(k) = lapdygreen(point,pts(k,:),ny(k,:))*ds(k);
            end
            u(i,j) = -mu'*muldygreen;
        end
    end
    
    for i = 1:2*M+1
        u(1,i)=boundaryFunc(x(1),y(i));
        u(2*M+1,i)=boundaryFunc(x(2*M+1),y(i));
        u(i,2*M+1)=boundaryFunc(x(i),y(2*M+1));
        u(i,2*M+2-ynumb(i))=boundaryFunc(x(i),y(2*M+2-ynumb(i)));
    end
    times(n) = toc();
    
    u=real(u);
    
    %计算误差
    e_2n = 0;ureal_2 = 0;
    for i = 1:2*M+1
        for j = 1:2*M+1
            if isnan(ureal(i,j))
            else
                e_2n = e_2n + (u(i,j)-ureal(i,j))^2;
                ureal_2 = ureal_2+ureal(i,j)^2;
            end
        end
    end
    e_2(n) = sqrt(e_2n/ureal_2);
    disp(e_2(n))
end
disp(times)

x = 200*(1:10); 
y = e_2;  

loglog(x, y, 'm*-');
grid on;
title('Error of the integral method with sin boundary');
xlabel('mesh number (log scale)');
ylabel('eror (log scale)');