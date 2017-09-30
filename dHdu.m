%% ����H�̓��͔���
% x        : [ x;dx ]       �i�ʒu�C���x�j
% u        : [ u;v;mu ]     �i����ʁC�_�~�[����ʁC���O�����W���搔�j
% lmd      : [ lmd1;lmd2 ]  �i�����ϐ�1�C�����ϐ�2�j
% sys      : a              �i�V�X�e���W���j
% sys      : b              �i�V�X�e���W���j
% cgmres   : [ r1;r2 ]      �i����ʂ̏d�݁C�_�~�[����ʂ̏d�݁j
% cgmres   : umax           �i����ʂ̍ő�l�j

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