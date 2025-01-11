e_2 = zeros(10,1);
times = zeros(10,1);
for n = 1:10   
    N = 100*n; % mesh number is 4*N
    
    % the boundary points
    pts = zeros(4*N,2); 
    
    pts(1:N,1) = ones(N,1);
    pts(N+1:2*N,1) = linspace(1-2/N,-1,N);
    pts(2*N+1:3*N,1) = -1*ones(N,1);
    pts(3*N+1:4*N,1) = linspace(-1+2/N,1,N);
    
    pts(1:N,2) = linspace(-1+2/N,1,N);
    pts(N+1:2*N,2) = ones(N,1);
    pts(2*N+1:3*N,2) = linspace(1-2/N,-1,N);
    pts(3*N+1:4*N,2) = -ones(N,1);
    
    % the normal derivitive
    ny = zeros(4*N,2);
    ny(1:N,1) = ones(N,1);
    ny(N+1:2*N,1) = zeros(N,1);
    ny(2*N+1:3*N,1) = -1*ones(N,1);
    ny(3*N+1:4*N,1) = zeros(N,1);
    
    ny(1:N,2) = zeros(N,1);
    ny(N+1:2*N,2) = ones(N,1);
    ny(2*N+1:3*N,2) = zeros(N,1);
    ny(3*N+1:4*N,2) = -ones(N,1);
    
    ny(N,:)=[1/2,1/2];
    ny(2*N,:)=[-1/2,1/2];
    ny(3*N,:)=[-1/2,-1/2];
    ny(4*N,:)=[1/2,-1/2];
    
    % boundary function
    boundaryFunc = @(x,y) x^2-y^2;
    
    % coefficients in the integral equation
    tic()
    alpha = (1/2)*ones(size(pts,1),1);
    for s = 1:4
        alpha(s*N)= 3/4;
    end
    dygreen = zeros(size(pts,1),size(pts,1));
    for i = 1:size(pts,1)
        for j = 1:size(pts,1)
            if i ~= j
                dygreen(i,j) = -lapdygreen(pts(i,:),pts(j,:),ny(j,:));
            else
                dygreen(i,j) = 0;
            end
        end
    end
    kernel = dygreen*8/size(pts,1);
    
    f = zeros(size(pts,1),1);
    for i = 1:size(pts,1)
        f(i) = boundaryFunc(pts(i,1),pts(i,2));
    end
    
    mu = solveieqn(kernel,alpha,f);
    
    % Integrate to get the solution
    M=10;
    x = linspace(-1,1,2*M+1);
    y = linspace(-1,1,2*M+1);
    u = zeros(2*M+1,2*M+1);
    muldygreen = zeros(size(pts,1),1);
    for i = 2:2*M
        for j = 2:2*M
            point = [x(i),y(j)];
            for k = 1:size(pts,1)
                muldygreen(k) = lapdygreen(point,pts(k,:),ny(k,:));
            end
            u(i,j) = -mu'*muldygreen*2/N;
        end
    end
    
    for i = 1:2*M+1
        u(1,i)=boundaryFunc(x(1),y(i));
        u(2*M+1,i)=boundaryFunc(x(2*M+1),y(i));
        u(i,1)=boundaryFunc(x(i),y(1));
        u(i,2*M+1)=boundaryFunc(x(i),y(2*M+1));
    end
    time = toc();

    %真解
    ureal = zeros(2*M+1,2*M+1);
    for i = 1:2*M+1
        for j = 1:2*M+1
            ureal(i,j) = boundaryFunc(x(i),y(j));
        end
    end
    
    %计算误差
    e_2(n) = norm(u-ureal)/norm(ureal);
    disp(e_2(n))
    times(n) = time;
end    

disp(e_2)

x = 100*(1:10); 
y = e_2;  

loglog(x, y, 'm*-'); 
grid on;
title('Error of the integral approach with square boundary');
xlabel('mesh number (log scale)');
ylabel('2-norm eror (log scale)');
