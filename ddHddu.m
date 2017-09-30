%% 函数Hの2階入力微分
% x        : [ x;dx ]       （位置，速度）
% u        : [ u;v;mu ]     （操作量，ダミー操作量，ラグランジュ乗数）
% lmd      : [ lmd1;lmd2 ]  （随伴変数1，随伴変数2）
% sys      : a              （システム係数）
% sys      : b              （システム係数）
% cgmres   : [ r1;r2 ]      （操作量の重み，ダミー操作量の重み）
% cgmres   : umax           （操作量の最大値）

function Huu = ddHddu( x, u, lmd, sys, cgmres )
    
    %{
    Huu = [ ...
        cgmres.r(1)+u(5)*2.0, 0, 0, 0, -cgmres.umax-cgmres.umin+u(1)*2.0, 0;
        0, cgmres.r(2)+u(6)*2.0, 0, 0, 0, -cgmres.umax-cgmres.umin+u(2)*2.0;
        0, 0, u(5)*2.0, 0, u(3)*2.0, 0;
        0, 0, 0, u(6)*2.0, 0, u(4)*2.0;
        -cgmres.umax-cgmres.umin+u(1)*2.0, 0, u(3)*2.0, 0, 0, 0;
        0, -cgmres.umax-cgmres.umin+u(2)*2.0, 0, u(4)*2.0, 0, 0;
    ];
    %}

    Huu = [ ...
        cgmres.r(1)+u(5)*2.0, 0, 0, 0, -cgmres.umax-cgmres.umin+u(1)*2.0, 0;
        0, cgmres.r(2)+u(6)*2.0, 0, 0, 0, -cgmres.umax-cgmres.umin+u(2)*2.0;
        0, 0, u(5)*2.0, 0, u(3)*2.0, 0;
        0, 0, 0, u(6)*2.0, 0, u(4)*2.0;
        -cgmres.umax-cgmres.umin+u(1)*2.0, 0, u(3)*2.0, 0, 0, 0;
        0, -cgmres.umax-cgmres.umin+u(2)*2.0, 0, u(4)*2.0, 0, 0;
    ];
end