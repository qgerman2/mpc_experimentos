clc;
clearvars;

A = [1.1 2; 0 0.95];
B = [0; 0.0787];
C = [-1 1];
D = 0;
Ts = 1;
sys = ss(A,B,C,D,Ts);
x0 = [0.5;-0.5]; % initial states at [0.5 -0.5]


t_unconstrained = 0:1:10;
u_unconstrained = zeros(size(t_unconstrained));

M = [A;A^2;A^3;A^4];
CONV = [B zeros(2,1) zeros(2,1) zeros(2,1);...
        A*B B zeros(2,1) zeros(2,1);...
        A^2*B A*B B zeros(2,1);...
        A^3*B A^2*B A*B B];

Q = diag([1, 1]);
R = 0.01;


Q_hat = blkdiag(Q,Q,Q,Q);
R_hat = blkdiag(R,R,R,R);
H = CONV'*Q_hat*CONV + R_hat;
F = CONV'*Q_hat*M;

K = H\F;
K_mpc = K(1,:);

figure(1);
Unconstrained_MPC = tf([-1 1])*feedback(ss(A,B,eye(2),0,Ts),K_mpc);
lsim(Unconstrained_MPC,'*',u_unconstrained,t_unconstrained,x0)
legend show

x(:,1) = x0;
u(:,1) = 0;
for i = 1:10
    x(:,i+1) = A*x(:,i)+B*-K_mpc*x(:,i);
    u(:,i+1) = -K_mpc*x(:,i);


end
hold on
plot(t_unconstrained, x')
figure(2);
plot(t_unconstrained, u);