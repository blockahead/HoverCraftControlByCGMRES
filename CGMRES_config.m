close all;
clear;


%% システム定義
M = 0.894;
Ig = 0.0125;
r = 0.0485;

sys.M = M;
sys.Ig = Ig;
sys.R = r;
sys.x0 = [-0.25;0.35;0;0;0;0;];

%% シミュレーション定義
tsim = 13; % シミュレーション時間 (s)


%% C/GMRESのコントローラ定義
dSamplingPeriod = 0.001;                % MATLAB Functionのサンプリング周期 (s)
cgmres.ht = 0.001;                      % 前進差分近似の時間幅 (s)
cgmres.zeta = 1000.0;                   % 操作量の安定化ゲイン (-)
cgmres.tf = 2.0;                        % 予測時間の最終値 (s)
cgmres.alpha = 1.0;                     % 予測時間の上昇速度ゲイン (-)
cgmres.dv = 5;                          % 予測時間の分割数 (-) （評価函数によって評価するポイントの数）

cgmres.x0 = [-0.25;0.35;0;0;0;0;];   % コントローラに与える初期状態

cgmres.q = [ 10;15;0.1;1;1;0.01; ];     % 状態に対する重み
cgmres.r = [ 1;1;0.01;0.01 ];           % 操作量に対する重み
cgmres.sf = cgmres.q;                   % 予測時間の最終状態に対する重み

cgmres.umin = -0.121;                   % 入力上限
cgmres.umax = 0.342;                    % 入力上限

desired = [ 0;0;0;0;0;0; ];

%% C/GMRESのコントローラ用計算
buff = c2d( ss( tf( [ cgmres.tf ], [ (1/cgmres.alpha), 1 ] ) ), dSamplingPeriod );
cgmres.T_outGain = buff.a;              % 予測時間差分方程式
cgmres.T_inGain = buff.b;               % 予測時間差分方程式
clearvars buff;

% 初期入力値の計算（Newton法）
lmd0 = dPhidx( cgmres.x0, desired, cgmres );
cgmres.u0 = [0;0;1;1;0.01;0.01;]; % Newton法の初期値

for cnt = 1:20
    cgmres.u0 = cgmres.u0 - ddHddu( cgmres.x0, cgmres.u0, lmd0, sys, cgmres ) \ dHdu( cgmres.x0, cgmres.u0, lmd0, sys, cgmres );
end

cgmres.len_x = length( cgmres.x0 );     % 状態の数
cgmres.len_u = length( cgmres.u0 );     % 操作量の数
cgmres.len_lmd = cgmres.len_x;          % 随伴変数の数


%% シミュレーションの実行と時間計測
tic;
sim( 'CGMRES' );
toc

%% 描画
hold on;
scatter( HoverOut(1,2), HoverOut(1,3) );
plot( HoverOut(:,2), HoverOut(:,3) );
len = 0.1;
direction = [ HoverOut(:,2) + len * cos( HoverOut(:,4) ), HoverOut(:,3) + len * sin( HoverOut(:,4) )];
for cnt = 1:20:length( HoverOut )
    plot( [ HoverOut(cnt,2);direction(cnt,1) ], [ HoverOut(cnt,3);direction(cnt,2) ] );
end

axis equal;

plot( [-0.3;0.3], [0;0] );
plot( [0;0], [-0.1;0.3] );

