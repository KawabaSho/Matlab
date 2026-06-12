clear,close all

% configuration
Iter  = 100;
dt    = 0.1;

% Initialization
YY    = zeros(1,Iter+1);
TT    = zeros(1,Iter+1);

YY(1,1) = 0;
% calculation
for i = 1 : Iter
    YY(1,i+1) = YY(1,i) + grad_sin(TT(1,i))*dt;
    TT(1,i+1) = dt*(i-1);
end
YY_Ref = sin(TT);

subplot(2,1,1)
plot(TT,YY,'-','Color','r')
hold on
plot(TT,YY_Ref,'-','Color','b')
hold off
legend(["Numerical Solution","Analytical Solution"])
xlabel("Time [s]")
ylabel("f(t) [-]")
subplot(2,1,2)
plot(TT,YY-YY_Ref,'-','Color','k')
xlabel("Time [s]")
ylabel("Error ( = Numerical - Analytical ) [-]")

% function
function dfdt = grad_sin(tt)
    dfdt = cos(tt);
end


