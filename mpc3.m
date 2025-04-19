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
% A = [A, B; zeros(m,n), eye(m)];
% B = [B; eye(m)];
% x0 = [x0; 0];
% 
% n = size(A,1);
% m = size(B,2);

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
q = diag([1, 1]);
r = 0;
Q = kron(eye(N), q);
R = kron(eye(N), r);

%% referencia
r = [0.1;0.1];
T = repmat(r, N, 1);

%% constraints
% Jan Maciejowski - Predictive control with constraints pag 82
% input
u_select = [1];
u_max = 4;
u_min = 4;
F = kron(tril(ones(N)), [u_select;-u_select]);
f = repmat([u_max; u_min], N, 1);

% input rate
du_select = [1];
du_max = 10;
du_min = 10;
A_du = kron(eye(N), [du_select;-du_select]);
b_du = repmat([du_max; du_min], N, 1);

% state
x_select = [
    1,0;
    0,1
];
x_max = [0.82;0.82];
x_min = [0.82;0.82];
Gamma = kron(eye(N), [x_select;-x_select]);
g = repmat([x_max;x_min], N, 1);

%% simular
% condicion inicial
x(:,1) = x0;
u = zeros(m, 1);
du = zeros(m, 1);

% optimizacion offline
H = Theta'*Q*Theta+R;
% H = (H+H')/2;

for k = 1:50
    % optimizacion online
    if (k > 1)
        epsilon = T - Psi*x(:,k) - Upsilon*u(:,k-1);
    else
        epsilon = T - Psi*x(:,k);
    end
    G = 2*Theta'*Q*epsilon;
    
    % constraints
    A_u = F;
    A_x = Gamma*Theta;
    if (k > 1)
        b_u = -F(:,1)*u(:,k-1)+f;
        b_x = -Gamma*(Psi*x(:,k)+Upsilon*u(:,k-1))+g;
    else
        b_u = f;
        b_x = -Gamma*Psi*x(:,k)+g;
    end

    % entrada mpc
    [opt, ~, exitflag] = quadprog(2*H, -G, [A_du; A_u; A_x],[b_du; b_u; b_x],[],[],[],[],[],options);
    if exitflag < 1
        disp("infeasible")
        break
    end
    du(:,k) = opt(1:m);
    if (k > 1)
        u(:,k) = u(:,k-1) + du(:,k);
    else
        u(:,k) = du(:,k);
    end

    % iterar
    x(:,k+1) = A*x(:,k)+B*u(:,k);
end

figure(1)
plot((x(1:2,:))')
figure(2);
plot(u);
figure(3);
plot(du)