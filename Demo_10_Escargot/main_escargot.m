%% ESCARGOT FOR YOU
% Algorithm: logarithmic spiral & frenet-Serret frame
% r,kn,kb,kr をいじってあなただけのESCARGOTを作りましょう

clear,close all

dt = 0.001;
Thf = 2*pi*4;
t = 0 : dt : Thf;

a = 1;
k = 0.1;
vz = 0.1;
xyz = escargot(t,a,k,vz);
r = 1-exp(-0 : 1/size(t,2) :1);
% r = 0 : 1/size(t,2) :1;
kn = 8;
kb = 1;
kr = 100;
ckr = cos(kr*t);
skr = sin(kr*t);
graxyz = grad_escargot(t,a,k,vz);
graxyz = graxyz./sqrt(sum(graxyz.^2));
Hessxyz = Hess_escargot(t,a,k,vz);
Hessxyz = Hessxyz./sqrt(sum(Hessxyz.^2));
for i = 1 : size(t,2)
    T = graxyz(1:3,i);
    P = xyz(1:3,i);
    N = Hessxyz(1:3,i);
    B = cross(T,N);
    xyz(1:3,i) = P + kn*r(i)*ckr(i)*N + kb*r(i)*skr(i)*B;
end

plot3(xyz(1,:),xyz(2,:),xyz(3,:),'Color',[0,0,0.7,0.1])
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.ZAxis.Visible = 'off';

function xyz = escargot(t,a,k,vz)
    xyz = [a*exp(k*t).*[cos(t);sin(t)]; -vz*t];
end
function xyz = grad_escargot(t,a,k,vz)
    xyz = [a*exp(k*t).*(k*[cos(t);sin(t)] + [-sin(t);cos(t)]); -vz.*ones(size(t))];
end
function xyz = Hess_escargot(t,a,k,~)
    xyz = [a*exp(k*t).*([k*k*cos(t)-sin(t)-k*sin(t)-cos(t);k*k*sin(t)+cos(t)+k*cos(t)-sin(t)]);
           zeros(size(t))];
end