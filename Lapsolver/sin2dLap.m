function [u] = sin2dLap(N,boundaryFunc)
    % Grid setup
    x = linspace(-1, 1, 2*N+1);
    y = linspace(-1, 1, 2*N+1);

    % 定义底边函数
    BottomBoundary = @(x) 1/2*(cos(x*pi)-1);
    %1/log(abs(x)/2)
    % 每列的格点个数
    ynumb = zeros(1,2*N+1,'int64');
    for i = 1:2*N+1
        if i ~= N+1
            ynumb(i) = 1+floor(N*(1-max(-1,BottomBoundary(x(i)))));
        else
            ynumb(i) = N+1;
        end
    end

    cumynumb = zeros(1,2*N+1,'int64'); % 前i-1列点的个数
    cumynumb(1) = 0;
    for i = 2:2*N+1
        cumynumb(i) = cumynumb(i-1)+ynumb(i-1);
    end
    inynumb = ynumb(1,2:2*N)-2*ones(1,2*N-1,'int64');
    totinnum = sum(inynumb);% 内部总格点数
    cuminynumb = zeros(1,2*N-1,'int64');
    cuminynumb(1) = 0;
    for i = 2:2*N-1
        cuminynumb(i) = cuminynumb(i-1)+inynumb(i-1);
    end
    
    % Initialize solution matrix
    u = zeros(2*N+1, 2*N+1);
    
    % nan part
    for i = 1:2*N+1
        for j = 1: 2*N-ynumb(i)+1
            u(i,j) = nan;
        end
    end
    
    % Construct the linear system
    % Interior grid points are indexed by 2:N+1 (skipping the boundary)
    A = zeros(totinnum,totinnum); 
    b = zeros(totinnum, 1);        % Right-hand side vector
    
    % Fill in A and b
    for i = 2:2*N
        for j = 2*N-ynumb(i)+3:2*N
            % Map 2D (i, j) index to 1D index
            k = cuminynumb(i-1)+j-2*N+ynumb(i)-2;
            
            % Coefficients for finite difference stencil
            A(k, k) = -4; % Center point

            if i == 2
                u(i-1,j) = boundaryFunc(x(i-1),y(j));
                b(k) = b(k) - u(i-1, j); % Left boundary
            elseif j < 2*N-ynumb(i-1)+3
                u(i-1,j) = boundaryFunc(x(i-1),y(j));
                b(k) = b(k) - u(i-1, j); % Left boundary
            else
                kleft = cuminynumb(i-2)+j-2*N+ynumb(i-1)-2;
                A(k, kleft) = 1; % Left neighbor
            end
            
            if i == 2*N
                u(i+1,j) = boundaryFunc(x(i+1),y(j));
                b(k) = b(k) - u(i+1, j); % Right boundary
            elseif j < 2*N-ynumb(i+1)+3
                u(i+1,j) = boundaryFunc(x(i+1),y(j));
                b(k) = b(k) - u(i+1, j); % Right boundary
            else
                kright = cuminynumb(i)+j-2*N+ynumb(i+1)-2;
                A(k,kright) = 1; % Right neighbor
            end
            
            if j == 2*N-ynumb(i)+3
                u(i,j-1) = boundaryFunc(x(i),y(j-1));
                b(k) = b(k) - u(i, j-1); % Bottom boundary
            else
                A(k, k-1) = 1;% Bottom neighbor
            end
            
            if j == 2*N
                u(i,j+1) = boundaryFunc(x(i),y(j+1));
                b(k) = b(k) - u(i, j+1); % Bottom boundary
            else
                A(k, k+1) = 1;% Bottom neighbor
            end
        end
    end
    
    % Solve the linear system
    uInterior = A \ b;
    
    % Populate the solution into the full matrix
    for i = 2:2*N
        for j = 2*N-ynumb(i)+3:2*N
            u(i, j) = uInterior(cuminynumb(i-1)+j-2*N+ynumb(i)-2);
        end
    end

end

