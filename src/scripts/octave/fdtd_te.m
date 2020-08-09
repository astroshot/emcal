% fdtd
% 二维FDTD  TE波圆柱仿真
% 定义常数
%-------------------------
c = 3.0E8;
mu = 1.2566E-6;
eps = 8.8542E-12;

f = 1E9;%频率
lambda = c/f;%波长
w_max = 2*pi*2*f;

% 定义FDTD网格
% dx\dy
ds = lambda/20;
%dt
dt = 0.5*ds/c;

dim_s = 5;
s_range = dim_s * lambda;
s_range = ceil(s_range / ds) * ds;%计算区域的长度
s_cells = s_range / ds;
cells = s_cells;%划分的网格数
nodes = cells+1;%采样点数

% 定义时间脉冲源
%----------------------------------------
dev_larger = 2;
dev = 1 / w_max * dev_larger;

dead = 4;
mean = dev * dead;

t = linspace(0, 9 * dev, 1000);
term = (t-mean);
pulse = (-1 / sqrt(2*pi) / dev^3) .* term;
pulse = pulse .* exp((-1 / 2 / dev^2) .* term .^2);
pulsenorm = max(pulse);

p = s_range / 2;
dead_s = 1.5;
dev_s = s_range / 4 / dead_s;
% 在空间划分采样点数
spacial = linspace(0,s_range,nodes);
taper = 1 / sqrt(2*pi) / dev_s * exp(-1 * (spacial - p) .^2 / 2 / dev_s^2);
taper = taper ./ max(taper);

figure(1);
plot(spacial,taper,'c');
title('脉冲源');
xlabel('x轴');
ylabel('Pluse');

% TE波的分量初始化
%----------------------------------------
Ex = zeros(nodes, nodes);
Ey = zeros(nodes, nodes);
Hz = zeros(nodes, nodes);

% 计算参数设置
%--------------------------------
done = 1;
n = 0;
F = 0;

c_mu = dt/mu/ds;
c_eps = dt/eps/ds;
c_MUR = (c*dt - ds)/(c*dt + ds);

frames=80;
figure(2);
axis([0 nodes 0 nodes 0 1]);
FDTD_M = moviein(frames);

while (done ~= 0),
    % 初始化
    Ex_o = Ex;
    Ey_o = Ey;
    Hz_o = Hz;

    % 电场值的中心差分
    Hterm = c_eps*(Hz(2:nodes-1,1:nodes) - Hz(1:nodes-2,1:nodes));
    Ex(2:nodes-1,1:nodes) = Ex(2:nodes-1,1:nodes) + Hterm;

    Hterm = c_eps*(Hz(1:nodes,1:nodes-2) - Hz(1:nodes,2:nodes-1));
    Ey(1:nodes,2:nodes-1) = Ey(1:nodes,2:nodes-1) + Hterm;

    %加入脉冲源
    t = n*dt;
    term = (t-mean);
    pulse = (-1/sqrt(2*pi)/dev^3)*term;
    pulse = pulse* exp((-1/2/dev^2)*term^2) / pulsenorm;
    source = pulse * ones(1,nodes) .* taper; % DBD correction

    Ex(2,1:nodes) = Ex(2,1:nodes) + source;

    % 4条边界线的处理
    % 左边
    Ex(1,1:nodes) = Ex_o(2,1:nodes) + c_MUR*(Ex(2,1:nodes) - Ex(1,1:nodes));
    % 上边
    Ey(1:nodes,nodes) = Ey_o(1:nodes,nodes-1) + c_MUR*(Ey(1:nodes,nodes-1) - Ey(1:nodes,nodes));
    % 右边
    Ex(nodes,1:nodes) = Ex_o(nodes-1,1:nodes) + c_MUR*(Ex(nodes-1,1:nodes) - Ex(nodes,1:nodes));
    % 下边
    Ey(1:nodes,1) = Ey_o(1:nodes,2) + c_MUR*(Ey(1:nodes,2) - Ey(1:nodes,1));

    % 磁场值的中心差分
    Eterm1 = c_mu*(Ex(2:nodes,1:nodes-1) - Ex(1:nodes-1,1:nodes-1));
    Eterm2 = c_mu*(Ey(1:nodes-1,1:nodes-1) - Ey(1:nodes-1,2:nodes));
    Hz(1:nodes-1,1:nodes-1) = Hz(1:nodes-1,1:nodes-1) + Eterm1 + Eterm2;

    % 描绘场分布
    cut = 3;
    figure(2);
    clf;
    if (abs(round(n/cut) - (n/cut)) == 0),  % 每5个时间步进行绘图
        F = F + 1;                          % 记录绘图的步数
        mesh(abs(Ex + i*Ey));               % 电场的幅值
        axis([0 nodes 0 nodes 0 1]);
        FDTD_M(:,F) = getframe;
    end;

    % 进行时间迭代，done=0时停止计算
    if (done ~= 0),
        n = n+1;                            % 记录迭代的时间步数
    end;
    if (n == cut*frames),
        done = 0;
    end;
end;

figure(3);
mesh(abs(Ex + i*Ey));
