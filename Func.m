%% èÛë‘ï˚íˆéÆ
% dx = f( x, u )
function dx = Func( x, u, sys )
    dx = [ ...
        x(4);
        x(5);
        x(6);
        (cos(x(3))*(u(1)+u(2)))/sys.M;
        (sin(x(3))*(u(1)+u(2)))/sys.M;
        (sys.R*(u(1)-u(2)))/sys.Ig;
    ];
      
end