clc;
clearvars;

options = optimoptions(@quadprog,'MaxIterations',1000,'ConstraintTolerance',1e-13,'Display','off');

%% parametros
n_BAT = 1;
Ts = 60; % s
C_E = 27000; % kJ
r_FC = 4.7E-8;
V_FC = 100;
N = 20;

%% sistema espacio estado discreto
% estado x      estado u
% SOC           P_FC
% H2
% P_M

A = [1, 0, -n_BAT/C_E*Ts; 
    0, 1, 0;
    0, 0, 1;
];
B = [
    n_BAT*Ts/C_E; 
    -r_FC*Ts/V_FC;
    0;
];

n = size(A,1);
m = size(B,2);

%% matrices propagacion estado
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
q = diag([10, 1, 0]);
r = 0;
Q = kron(eye(N), q);
R = kron(eye(N), r);

%% referencia
r = [0.9;0.9;0];
T = repmat(r, N, 1);

%% constraints
% input
u_select = [1];
u_max = 128;
u_min = 0;
F = kron(tril(ones(N)), [u_select;-u_select]);
f = repmat([u_max; u_min], N, 1);

% input rate
du_select = [1];
du_max = 15*Ts;
du_min = 15*Ts;
A_du = kron(eye(N), [du_select;-du_select]);
b_du = repmat([du_max; du_min], N, 1);

% state
x_select = [
    1,0,0;
    0,1,0
];
x_max = [1;1];
x_min = [0;0];
Gamma = kron(eye(N), [x_select;-x_select]);
g = repmat([x_max;x_min], N, 1);

%% simular
% condicion inicial
x(:,1) = [0.5, 0.5, 100];
u = zeros(m, 1);
du = zeros(m, 1);

% optimizacion offline
H = Theta'*Q*Theta+R;
H = (H+H')/2;

for k = 1:20
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
        b_u = -F(:,1) * u(:,k-1) + f;
        b_x = -Gamma * (Psi*x(:,k) + Upsilon*u(:,k-1)) + g;
    else
        b_u = f;
        b_x = -Gamma * Psi*x(:,k) + g;
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
ylabel("x")
figure(2);
plot(u);
ylabel("u")
figure(3);
plot(du)
ylabel("du")