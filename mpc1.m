clc;
clearvars;

options = optimoptions(@quadprog,'MaxIterations',1000,'ConstraintTolerance',1e-13,'Display','off');

A = [1.1 2; 0 0.95];
B = [0; 0.0787];
C = [-1 1];
D = 0;
Ts = 1;

x0 = [0.5;-0.5];

N = 10;
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
% F = zeros(n*N,n);
% G = zeros(n*N,m*(N-1));
% 
% % form row by row
% for i = 1:N
%     % F
%     F(n*(i-1)+(1:n),:) = A^i;
% 
%     % G
%     % form element by element
%     for j = 1:i
%         G(n*(i-1)+(1:n),m*(j-1)+(1:m)) = A^(i-j)*B;
%     end
% end

%% costos
q = diag([1, 1]);
r = 0;
Q = kron(eye(N), q);
R = kron(eye(N), r);

%% referencia
r = [0;0];
T = repmat(r, N, 1);

%% simular
% condicion inicial
x(:,1) = x0;
u = zeros(m, 1);
p_u = zeros(m, 1);

% optimizacion offline
% H = 2*G'*Q*G + 2*R;
% L = 2*G'*Q*F;
H = Theta'*Q*Theta+R;

figure(3)
for k = 1:10
    epsilon = T - Psi*x(:,k) - Upsilon*p_u;
    G = 2*Theta'*Q*epsilon;
    
    %opt = 1/2*inv(H)*G;
    [opt, ~, exitflag] = quadprog(2*H, -G, [],[],[],[],[],[],[],options);
    hold on
    plot(k:k+N-1, opt)

    u(:,k) = p_u + opt(1);

    x(:,k+1) = A*x(:,k)+B*u(:,k);
    p_u = u(:,k);
end

figure(1)
plot(x')
figure(2);
plot(u)