clear, close all

N = 1000;
tf = 100;
dt = tf / N;

tt = zeros(1,N+1);
X = zeros(2,N+1);

% initialization
X(:,1) = [0;0];
u = 1; % 1 [N]

m = 3;
c = 2;
k = 1;

ome = sqrt(k/m);
zeta = c/m/2/ome;

A = [0, 1; -ome^2, -2*zeta*ome;];
B = [0; 1];

func = @(x,u)grad_func(x,u, A, B);

for i = 1 : N
    tt(i+1) = dt*i;
    X(:,i+1) = X(:, i) + func(X(:, i), u)*dt;
    pause;
end

% Draw
plot(tt, X(1,:))
hold on
plot(tt, X(2,:))
hold off
legend(["x","dot x"])




function dxdt = grad_func(x_k, u, A, B)
    dxdt = A*x_k + B*u;
end


