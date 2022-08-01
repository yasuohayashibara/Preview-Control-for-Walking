% preview control for walking

clear;
											%足の着地時間(s)，x方向の位置(m)，y方向の位置(m)
foot = [0 0 0; 0.6 0.1 0.06; 0.9 0.2 -0.06; 1.2 0.3 0.06; 1.5 0.4 -0.06; 1.8 0.5 0.06; 2.4 0.6 -0.06; 3.0 0.7 0.0; 100 0 0];
forward_period = 1.0;						% 予見制御の時間(s)
calculate_period = 4.0;						% 歩行パターンを生成する期間(s)
dt = 0.01;									% サンプリングタイム (s)
zh = 0.27;									% 重心位置 (m)
g  = 9.8;									% 重力加速度 (m/s^2)
A = [0 1 0; 0 0 1; 0 0 0];
B = [0; 0; 1];
C = [1 0 -zh/g];
D = 0;
sys = ss(A, B, C, D);
sys_d = c2d(sys, dt);
% [A_d, B_d, C_d, D_d] = sys2ss(sys_d);
[A_d, B_d, C_d, D_d] = ssdata(sys_d);
E_d = [dt; 1; 0];

Zero = [0; 0; 0];
Phai = [1 -C_d*A_d; Zero A_d];
G = [-C_d*B_d; B_d];
GR = [1; Zero];
Gd = [-C_d*E_d; E_d];

Q = zeros(4);
Q(1) = 10^8;								% 最適制御の重み係数　大：ZMPが理想に近づく
H = 1;										% 同じく重み係数　大：入力が少なくなる．

P = dare(Phai, G, Q, H);
F = -(H+G'*P*G)^(-1)*G'*P*Phai;

x = [0; 0; 0];
y = [0; 0; 0];
xp = x;
yp = x;

t = 0:dt:calculate_period;

i = 1;										% 目標ZMPの生成
n = 1;
for tt = 0:dt:calculate_period+forward_period+1
	if (tt == foot(n,1))
		prefx(i) = foot(n,2);
		prefy(i) = foot(n,3);
		n = n + 1;
	else
		prefx(i) = prefx(i-1);
		prefy(i) = prefy(i-1);
    end
	i = i + 1;
end

i = 0;
ux = 0;
uy = 0;

xi = (eye(4)-G*(H+G'*P*G)^(-1)*G'*P)*Phai;

for tt = t
	i = i + 1;
%	fd = (H+G'*P*G)^(-1)*G'*(xi')^1*P*Gd;	% 外乱除去用　未使用
	px = C_d*x;
	py = C_d*y;
	ex = prefx(i) - px;
	ey = prefy(i) - py;
	X = [ex; x - xp];
	Y = [ey; y - yp];
	xp = x;
	yp = y;
	dux = F * X;
	j = 0;
	for ttt = tt : dt : (tt + forward_period)
		j = j + 1;
		if (prefx(i+j) - prefx(i+j-1)) ~= 0
			f  = -(H+G'*P*G)^(-1)*G'*(xi')^(j-1)*P*GR;
			dux = dux + f * (prefx(i+j) - prefx(i+j-1));
        end
    end
	ux = ux + dux;
	duy = F * Y;
	j = 0;
	for ttt = tt : dt : (tt + forward_period)
		j = j + 1;
		if (prefy(i+j) - prefy(i+j-1)) ~= 0
			f  = -(H+G'*P*G)^(-1)*G'*(xi')^(j-1)*P*GR;
			duy = duy + f * (prefy(i+j) - prefy(i+j-1));
        end
    end
	uy = uy + duy;

	dx = 0;
	dy = 0;
	x = A_d * x + B_d * ux + E_d * dx * dt;
	y = A_d * y + B_d * uy + E_d * dy * dt;
	x0(i) = x(1);
	y0(i) = y(1);
	x1(i) = prefx(i);
	y1(i) = prefy(i);
	x2(i) = px;
	y2(i) = py;
end

subplot(2,1,1);
plot(x0, y0, "o", x1, y1, x2, y2);
subplot(2,1,2);
plot(t, x0, t, x1, t, x2, t, y0, t, y1, t, y2);
