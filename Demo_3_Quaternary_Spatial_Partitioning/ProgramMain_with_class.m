clear,close all
% Quaternary Tree 2024-03-30
% Treeで空間レベルとその空間内での番号が分かります
% 空間の可視化はまた後日



N = 10;
ParticleSize = 0.015;
X = rand(2,N)*1/(1+ParticleSize*2)+ParticleSize;
% X = [0.370282112094186,	0.699088357087138;
% 0.764508608875511,	0.671364320606215];
% X = [0.81;0.7];
% load X.mat
xy_range = [0,1];
dX = xy_range(2)-xy_range(1);


PartitionLevel = 5;

QSP = class_Quaternary_Spatial_Partitiong(PartitionLevel,xy_range,ParticleSize);
tic
for k = 1 : N
    QSP.AddObject(X(1,k),X(2,k))
end
toc
% [level, num]



% graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
DrawCircle(X(1,1),X(2,1),ParticleSize,1)
hold on
for i = 1 : N
    DrawCircle(X(1,i),X(2,i),ParticleSize,i+1)
    text(X(1,i),X(2,i),num2str(i))
end
U = (xy_range(2)-xy_range(1))/2^PartitionLevel;
DrawLine(1,U)
DrawBox(0.5,0.5,1)

% Tree
for i = 1 : N
    LId = double(QSP.Object2Tree{i});
    LL = dX*0.5^LId(1);
    [x,y] = Box2Position(LId(2),LL);
    DrawBox(x,y,LL)
end

axis([-0.05,1+0.05,-0.05,1+0.05])
axis square 







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DrawCircle(x,y,r,n)
    t = linspace(0,2*pi,100);
    plot(r*sin(t)+x,r*cos(t)+y,'Color',getColor(n))
end
function DrawBox(x,y,L)
    s = [-L,-L,L,L,-L; -L,L,L,-L,-L;]*0.5 + [x;y];
    plot(s(1,:),s(2,:),'k')
end
function DrawLine(maxL,U)
    num = maxL/U;
    for i = 1:num+1
        k = i - 1;
        d = 0.6;
        if k*U == 0.5
            d = 2*d;
        end
        xline(k*U,'Color',[0.9,0.9,0.9],'LineWidth',d);
        yline(k*U,'Color',[0.9,0.9,0.9],'LineWidth',d);
    end
end
function [x,y] = Box2Position(Num,dX)
    x = double(BitGet32(Num));
    y = double(BitGet32(bitshift(Num,-1)));
    x = (x+0.5)*dX;
    y = (y+0.5)*dX;
end
function intR = BitGet32(n)
    % 16bit
    intR = uint16(0);
    for i = 1 : 16
        intR = bitset(intR,i,bitget(n,2*i-1));
    end
end


function s = getColor(n)
    n = mod(n,8)+1;
    p =[0.85,0.33,0.10;
        0.00,0.45,0.74;
        0.47,0.67,0.19;
        0.85,0.33,0.10;
        0.49,0.18,0.56;
        0.64,0.08,0.18;
        0.93,0.69,0.13;
        0.49,0.18,0.56;];
    s = p(n,:);
end