%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load control
clear;
dist = [0.6, 0.05, 3.14/2];
max_step_x = 0.1;
max_step_y = 0.05;
max_step_w = 0.2;
period = 0.3;
foot_y = 0.06;
step_num_x = fix(abs(dist(1))/max_step_x);
step_num_y = fix(abs(dist(2))/max_step_y);
step_num_w = fix(abs(dist(3))/max_step_w);
step_num = max([step_num_x, step_num_y, step_num_w]);
if (step_num_x == step_num)
	step_x = sign(dist(1))*min(abs(dist(1)), max_step_x);
	step_y = dist(2)/abs(dist(1))*max_step_x;
	step_w = dist(3)/abs(dist(1))*max_step_x;
elseif (step_num_y == step_num)
	step_x = dist(1)/abs(dist(2))*max_step_y;
	step_y = sign(dist(2))*min(abs(dist(2)), max_step_y);
	step_w = dist(3)/abs(dist(2))*max_step_y;
elseif (step_num_w == step_num)
	step_x = dist(1)/abs(dist(3))*max_step_w;
	step_y = dist(2)/abs(dist(3))*max_step_w;
	step_w = sign(dist(3))*min(abs(dist(3)), max_step_w);
endif

foot(1,:) = [0 0 0];
foot(2,:) = [period 0 foot_y];
rot = 0;
for i = 2:step_num+1
	shift_y = foot_y * (rem(i,2)*2-1)*2;
	foot_rd = [step_x, step_y+shift_y];
	foot_fd = [foot_rd(1)*cos(rot)-foot_rd(2)*sin(rot), foot_rd(1)*sin(rot)+foot_rd(2)*cos(rot)];
	rot = rot + step_w
	foot(i+1,:) = foot(i,:)+[period, foot_fd(1), foot_fd(2)];
endfor
foot_rd = [dist(1)-step_x*step_num, dist(2)-step_y*step_num-shift_y];
foot_fd = [foot_rd(1)*cos(rot)-foot_rd(2)*sin(rot), foot_rd(1)*sin(rot)+foot_rd(2)*cos(rot)];
foot(step_num+3,:) = foot(step_num+2,:) + [period, foot_fd(1), foot_fd(2)];
rot = dist(3)
foot_rd = [0,shift_y/2];
foot_fd = [foot_rd(1)*cos(rot)-foot_rd(2)*sin(rot), foot_rd(1)*sin(rot)+foot_rd(2)*cos(rot)];
foot(step_num+4,:) = foot(step_num+3,:) + [period, foot_fd(1), foot_fd(2)];
foot(step_num+5,:) = [100, 0, 0];

forward_period = 1.0;						% 予見制御の時間(s)
calculate_period = 4.0;						% 歩行パターンを生成する期間(s)
dt = 0.01;									% サンプリングタイム (s)
zh = 0.27;									% 重心位置 (m)
g  = 9.8;									% 重力加速度 (m/s^2)
A = [0 1 0; 0 0 1; 0 0 0];					% 式(5.8)
B = [0; 0; 1];
C = [1 0 -zh/g];							% 式(5.9)
D = 0;
sys = ss(A, B, C, D);
sys_d = c2d(sys, dt);						% 離散化
if (substr(version,1,3)=="3.2")
	[A_d, B_d, C_d, D_d] = sys2ss(sys_d);	% 式(5.10)
else
	[A_d, B_d, C_d, D_d] = ssdata(sys_d);
endif
E_d = [dt; 1; 0];							% 式(5.11)

Zero = [0; 0; 0];							% 式(5.12),(5.13)
Phai = [1 -C_d*A_d; Zero A_d];
G = [-C_d*B_d; B_d];
GR = [1; Zero];
Gd = [-C_d*E_d; E_d];

Q = zeros(4);
Q(1) = 10^8;								% ★★★最適制御の重み係数　大：ZMPが理想に近づく
H = 1;										% 同じく重み係数　大：入力が少なくなる．

P = dare(Phai, G, Q, H);					% 離散時間系でのリカッチ方程式の解
F = -(H+G'*P*G)^(-1)*G'*P*Phai;				% ★★★フィードバックゲイン 式(2.13)

x = [0; 0; 0];
y = [0; 0; 0];
xp = x;
yp = x;

t = 0:dt:calculate_period;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 1;										% 目標ZMP波形の生成
n = 1;										% 足の着地時間・位置 foot → 目標ZMPパターン prefx, prefy
for tt = 0:dt:calculate_period+forward_period+1
	if (abs(tt - foot(n,1))<(dt/2))
		prefx(i) = foot(n,2);
		prefy(i) = foot(n,3);
		n = n + 1;
	else
		prefx(i) = prefx(i-1);
		prefy(i) = prefy(i-1);
	endif
	i = i + 1;
endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 0;
ux = uy = 0;

xi = (eye(4)-G*(H+G'*P*G)^(-1)*G'*P)*Phai;

for tt = t
	i = i + 1;
%	fd = (H+G'*P*G)^(-1)*G'*(xi')^1*P*Gd;	% 外乱除去用　未使用
	px = C_d*x;								% 実際のZMP位置
	py = C_d*y;
	ex = prefx(i) - px;						% 目標ZMP位置と実際との差
	ey = prefy(i) - py;
	X = [ex; x - xp];						% 式(2.5),(2,6)
	Y = [ey; y - yp];
	xp = x;									% 差分を計算するため，一つ前の値を保存
	yp = y;
	dux = F * X;							% 式(2.8)
	j = 0;
	for ttt = tt : dt : (tt + forward_period)
		j = j + 1;
		if (prefx(i+j) - prefx(i+j-1)) != 0
			f  = -(H+G'*P*G)^(-1)*G'*(xi')^(j-1)*P*GR;		% ★★★式(2.14)
			dux = dux + f * (prefx(i+j) - prefx(i+j-1));	% 式(2.12)
		endif
	endfor
	ux += dux;
	duy = F * Y;							% y方向の予見制御
	j = 0;
	for ttt = tt : dt : (tt + forward_period)
		j = j + 1;
		if (prefy(i+j) - prefy(i+j-1)) != 0
			f  = -(H+G'*P*G)^(-1)*G'*(xi')^(j-1)*P*GR;		% ★★★式(2.14)
			duy = duy + f * (prefy(i+j) - prefy(i+j-1));	% 式(2.12)
		endif
	endfor
	uy += duy;

	dx = 0;									% 外乱（今回は無し）
	dy = 0;
	x = A_d * x + B_d * ux + E_d * dx * dt;	% 重心位置の計算
	y = A_d * y + B_d * uy + E_d * dy * dt;
	x0(i) = x(1);							% 重心位置の描画用
	y0(i) = y(1);
	x1(i) = prefx(i);						% 目標ZMP位置の描画用
	y1(i) = prefy(i);
	x2(i) = px;								% 実際のZMP位置の描画用
	y2(i) = py;
endfor
subplot(2,1,1);								% グラフの描画
plot(x0, y0, "o", x1, y1, x2, y2);
subplot(2,1,2);
plot(t, x0, t, x1, t, x2, t, y0, t, y1, t, y2);
