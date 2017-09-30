%% 函数Hの状態微分
% x        : [ x;dx ]       （位置，速度）
% lmd      : [ lmd1;lmd2 ]  （随伴変数1，随伴変数2）
% sys      : a              （システム係数）
% sys      : b              （システム係数）
% cgmres   : [ q1;q2 ]      （位置の重み，速度の重み）

function Hx = dHdx( x, u, lmd, d, sys, cgmres )

    %{
    Hx = [ ...
        cgmres.q(1)*x(1);
        cgmres.q(2)*x(2);
        cgmres.q(3)*x(3)+(cos(x(3))*(u(1)+u(2))*lmd(5))/sys.M-(sin(x(3))*(u(1)+u(2))*lmd(4))/sys.M;
        cgmres.q(4)*x(4)+lmd(1);
        cgmres.q(5)*x(5)+lmd(2);
        cgmres.q(6)*x(6)+lmd(3);
    ];
    %}

    Hx = [ ...
        (d(1)*2.0-x(1)*2.0)*cgmres.q(1)*(-1.0/2.0);
        (d(2)*2.0-x(2)*2.0)*cgmres.q(2)*(-1.0/2.0);
        (d(3)*2.0-x(3)*2.0)*cgmres.q(3)*(-1.0/2.0)+(cos(x(3))*(u(1)+u(2))*lmd(5))/sys.M-(sin(x(3))*(u(1)+u(2))*lmd(4))/sys.M;
        lmd(1)-(d(4)*2.0-x(4)*2.0)*cgmres.q(4)*(1.0/2.0);
        lmd(2)-(d(5)*2.0-x(5)*2.0)*cgmres.q(5)*(1.0/2.0);
        lmd(3)-(d(6)*2.0-x(6)*2.0)*cgmres.q(6)*(1.0/2.0);
    ];
end