e_2 = zeros(10,1);
times = zeros(10,1);
for k = 1:10 % 网格数目
    N = 10*k;
    boundaryFunc = @(x, y) sin(x)*cosh(y);
    tic()
    u = SquareLap(N, boundaryFunc);
    time = toc();
    fprintf('4-point scheme time: %10.4e (s) \n',time)
    x = linspace(-1, 1, N+1);
    y = linspace(-1, 1, N+1);
    ureal = zeros(N+1,N+1);
    for i = 1:N+1
        for j = 1:N+1
            ureal(i,j) = boundaryFunc(x(i),y(j));
        end
    end
    
    % 计算误差
    e_2k = norm(u-ureal)/norm(ureal);
    
    e_2(k) = e_2k;
    times(k) = time;
end
    
disp(e_2)

x = 10*(1:10); 
y = e_2;  

loglog(x, y, 'm*-'); 
grid on;
title('Error of the 4-points method with square boundary');
xlabel('mesh number (log scale)');
ylabel('2-norm eror (log scale)');

