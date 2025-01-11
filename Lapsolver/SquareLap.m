function [u] = SquareLap(N,boundaryFunc)
% solveLaplaceDirichlet - Solve Laplace's equation with Dirichlet boundary
    % conditions on the domain [-1, 1] x [-1, 1].
    %
    % Inputs:
    %   N - Number of interior grid points along one dimension (total grid will be (N+2)x(N+2))
    %   boundaryFunc - Function handle for boundary conditions: boundaryFunc(x, y)
    %
    % Outputs:
    %   u - Solution matrix (including boundary values)
    
    % Grid setup
    x = linspace(-1, 1, N+1);
    y = linspace(-1, 1, N+1);
    
    % Initialize solution matrix
    u = zeros(N+1, N+1);
    
    % Apply boundary conditions
    for i = 1:N+1
        u(i, 1) = boundaryFunc(x(i), y(1));       % Bottom boundary
        u(i, N+1) = boundaryFunc(x(i), y(N+1));   % Top boundary
        u(1, i) = boundaryFunc(x(1), y(i));       % Left boundary
        u(N+1, i) = boundaryFunc(x(N+1), y(i));   % Right boundary
    end
    
    % Construct the linear system
    % Interior grid points are indexed by 2:N+1 (skipping the boundary)
    A = zeros((N-1)^2, (N-1)^2); % Sparse matrix for efficiency
    b = zeros((N-1)^2, 1);        % Right-hand side vector
    
    % Fill in A and b
    for j = 1:N-1
        for i = 1:N-1
            % Map 2D (i, j) index to 1D index
            k = (j-1)*(N-1) + i;
            
            % Coefficients for finite difference stencil
            A(k, k) = -4; % Center point
            if i > 1
                A(k, k-1) = 1; % Left neighbor
            else
                b(k) = b(k) - u(1, j+1); % Left boundary
            end
            
            if i < N-1
                A(k, k+1) = 1; % Right neighbor
            else
                b(k) = b(k) - u(N+1, j+1); % Right boundary
            end
            
            if j > 1
                A(k, k-N+1) = 1; % Bottom neighbor
            else
                b(k) = b(k) - u(i+1, 1); % Bottom boundary
            end
            
            if j < N-1
                A(k, k+N-1) = 1; % Top neighbor
            else
                b(k) = b(k) - u(i+1, N+1); % Top boundary
            end
        end
    end
    
    % Solve the linear system
    uInterior = A \ b;
    
    % Populate the solution into the full matrix
    for j = 1:N-1
        for i = 1:N-1
            u(i+1, j+1) = uInterior((j-1)*(N-1) + i);
        end
    end

end

