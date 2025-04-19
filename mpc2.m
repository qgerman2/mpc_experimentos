clc;
clearvars;

options = optimoptions(@quadprog,'MaxIterations',1000,'ConstraintTolerance',1e-13,'Display','off');

A = [1.1 2; 0 0.95];
B = [0; 0.0787];
x0 = [0.5;-0.5];

N = 10;
n = size(A,1);
m = size(B,2);

% extender estado con input previo
A = [A, B; zeros(m,n), eye(m)];
B = [B; eye(m)];
x0 = [x0; 0];

n = size(A,1);
m = size(B,2);

%% matrices propagacion estado
% Jan Maciejowski - Predictive control with constraints pag 55
Psi = zeros(N*n,n);
for i = 1:N
    Psi(n*(i-1)+(1:n),:) = A^i;
end
Upsilon = zeros(N*n,m);
for i = 1:N
    for j = 1:i
        Upsilon(n*(i-1)+(1:n),:) = Upsilon(n*(i-1)+(1:n),:) + ...
            A^(j-1) * B;
    end
end
Theta = zeros(N*n,N*m);
for i = 1:N
    for j = 1:i
        for k = 1:i-j+1
            Theta(n*(i-1)+(1:n),m*(j-1)+(1:m)) = Theta(n*(i-1)+(1:n),m*(j-1)+(1:m)) + ...
                A^(k-1) * B;
        end
    end
end

%% costos
q = diag([1, 1, 0]);
r = 0;
Q = kron(eye(N), q);
R = kron(eye(N), r);

%% referencia
r = [0;0;0];
T = repmat(r, N, 1);

%% simular
% condicion inicial
x(:,1) = x0;
u = zeros(m, 1);

% optimizacion offline
H = Theta'*Q*Theta+R;

for k = 1:10
    % optimizacion online
    epsilon = T - Psi*x(:,k) - Upsilon*x(n,k);
    G = 2*Theta'*Q*epsilon;
    
    [opt, ~, exitflag] = quadprog(2*H, -G, [],[],[],[],[],[],[],options);

    u(:,k) = x(n,k) + opt(1);
    x(:,k+1) = A*x(:,k)+B*u(:,k);
end

figure(1)
plot((x(1:2,:))')
figure(2);
plot(u)