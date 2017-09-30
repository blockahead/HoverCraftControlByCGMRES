%% 函数Hの入力微分
% x        : [ x;dx ]       （位置，速度）
% u        : [ u;v;mu ]     （操作量，ダミー操作量，ラグランジュ乗数）
% lmd      : [ lmd1;lmd2 ]  （随伴変数1，随伴変数2）
% sys      : a              （システム係数）
% sys      : b              （システム係数）
% cgmres   : [ r1;r2 ]      （操作量の重み，ダミー操作量の重み）
% cgmres   : umax           （操作量の最大値）

function Hu = dHdu( x, u, lmd, sys, cgmres )
    
    %{
    Hu = [ ...
        cgmres.r(1)*u(1)-u(5)*(cgmres.umax+cgmres.umin-u(1)*2.0)+(cos(x(3))*lmd(4))/sys.M+(sin(x(3))*lmd(5))/sys.M+(sys.R*lmd(6))/sys.Ig;
        cgmres.r(2)*u(2)-u(6)*(cgmres.umax+cgmres.umin-u(2)*2.0)+(cos(x(3))*lmd(4))/sys.M+(sin(x(3))*lmd(5))/sys.M-(sys.R*lmd(6))/sys.Ig;
        u(3)*u(5)*2.0-cgmres.r(3);
        u(4)*u(6)*2.0-cgmres.r(4);
        u(3)^2-(cgmres.umax*(1.0/2.0)-cgmres.umin*(1.0/2.0))^2+(cgmres.umax*(1.0/2.0)+cgmres.umin*(1.0/2.0)-u(1))^2;
        u(4)^2-(cgmres.umax*(1.0/2.0)-cgmres.umin*(1.0/2.0))^2+(cgmres.umax*(1.0/2.0)+cgmres.umin*(1.0/2.0)-u(2))^2;
    ];
    %}

    Hu = [ ...
        cgmres.r(1)*u(1)-u(5)*(cgmres.umax+cgmres.umin-u(1)*2.0)+(cos(x(3))*lmd(4))/sys.M+(sin(x(3))*lmd(5))/sys.M+(sys.R*lmd(6))/sys.Ig;
        cgmres.r(2)*u(2)-u(6)*(cgmres.umax+cgmres.umin-u(2)*2.0)+(cos(x(3))*lmd(4))/sys.M+(sin(x(3))*lmd(5))/sys.M-(sys.R*lmd(6))/sys.Ig;
        u(3)*u(5)*2.0-cgmres.r(3);
        u(4)*u(6)*2.0-cgmres.r(4);
        u(3)^2-(cgmres.umax*(1.0/2.0)-cgmres.umin*(1.0/2.0))^2+(cgmres.umax*(1.0/2.0)+cgmres.umin*(1.0/2.0)-u(1))^2;
        u(4)^2-(cgmres.umax*(1.0/2.0)-cgmres.umin*(1.0/2.0))^2+(cgmres.umax*(1.0/2.0)+cgmres.umin*(1.0/2.0)-u(2))^2;
    ];

end