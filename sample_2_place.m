clear, close all
%%
N = 1000;
tf = 12;
dt = tf / N;

tt = zeros(1,N+1);
X = zeros(2,N+1);

% initialization
X(:,1) = [0;0];


m = 1;
c = 0.5;
k = 1;

ome = sqrt(2);
zeta = 1/2/ome;

A = [0, 1; -ome^2, -2*zeta*ome;];
B = [0; 2];
u_c = 1;

% place 
p = [-1-2*i,-1+2*i];
K = place(A,B,p);


func = @(x,u)grad_func(x,u, A, B);

for i = 1 : N
    
    x = X(:, i);

    % u =  - K*x;
    u = u_c - K*x;
    % u = u_c;

    k1 = func(x,u)*dt;
    k2 = func(x+0.5*k1,u)*dt;
    k3 = func(x+0.5*k2,u)*dt;
    k4 = func(x+k3,u)*dt;

    X(:,i+1) = x + (k1+(k2+k3)*2+k3)/6;
    tt(i+1) = dt*i;

    % pause;
end
%%
% Draw
figure("Name","sample1")
subplot(2,1,1)
hold on 
plot(tt, X(1,:))
plot([0,12],[1,1],'--','Color',[.1,.1,.1])
ylabel('x [m]')
hold off
subplot(2,1,2)
hold on 
plot(tt, X(2,:))
plot([0,12],[0,0],'--','Color',[.1,.1,.1])
hold off
ylabel('v [m/s]')
xlabel('time [s]')


function dxdt = grad_func(x_k, u, A, B)
    dxdt = A*x_k + B*u;
end


