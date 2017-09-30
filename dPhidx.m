%% ヘッセ行列の状態微分
% x        : [ x;dx ]       （位置，速度）
% cgmres   : [ sf1;sf2 ]    （時刻Tでの位置の重み，時刻Tでの速度の重み）

function Phix = dPhidx( x, d, cgmres )

    %{
    Phix = [ ...
        cgmres.sf(1)*x(1);
        cgmres.sf(2)*x(2);
        cgmres.sf(3)*x(3);
        cgmres.sf(4)*x(4);
        cgmres.sf(5)*x(5);
        cgmres.sf(6)*x(6);
    ];
    %}
    
    Phix = [ ...
        (d(1)*2.0-x(1)*2.0)*cgmres.sf(1)*(-1.0/2.0);
        (d(2)*2.0-x(2)*2.0)*cgmres.sf(2)*(-1.0/2.0);
        (d(3)*2.0-x(3)*2.0)*cgmres.sf(3)*(-1.0/2.0);
        (d(4)*2.0-x(4)*2.0)*cgmres.sf(4)*(-1.0/2.0);
        (d(5)*2.0-x(5)*2.0)*cgmres.sf(5)*(-1.0/2.0);
        (d(6)*2.0-x(6)*2.0)*cgmres.sf(6)*(-1.0/2.0);
    ];

end