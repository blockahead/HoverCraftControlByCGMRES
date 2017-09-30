close all;
clear;


%% �V�X�e����`
M = 0.894;
Ig = 0.0125;
r = 0.0485;

sys.M = M;
sys.Ig = Ig;
sys.R = r;
sys.x0 = [-0.25;0.35;0;0;0;0;];

%% �V�~�����[�V������`
tsim = 13; % �V�~�����[�V�������� (s)


%% C/GMRES�̃R���g���[����`
dSamplingPeriod = 0.001;                % MATLAB Function�̃T���v�����O���� (s)
cgmres.ht = 0.001;                      % �O�i�����ߎ��̎��ԕ� (s)
cgmres.zeta = 1000.0;                   % ����ʂ̈��艻�Q�C�� (-)
cgmres.tf = 2.0;                        % �\�����Ԃ̍ŏI�l (s)
cgmres.alpha = 1.0;                     % �\�����Ԃ̏㏸���x�Q�C�� (-)
cgmres.dv = 5;                          % �\�����Ԃ̕����� (-) �i�]�������ɂ���ĕ]������|�C���g�̐��j

cgmres.x0 = [-0.25;0.35;0;0;0;0;];   % �R���g���[���ɗ^���鏉�����

cgmres.q = [ 10;15;0.1;1;1;0.01; ];     % ��Ԃɑ΂���d��
cgmres.r = [ 1;1;0.01;0.01 ];           % ����ʂɑ΂���d��
cgmres.sf = cgmres.q;                   % �\�����Ԃ̍ŏI��Ԃɑ΂���d��

cgmres.umin = -0.121;                   % ���͏��
cgmres.umax = 0.342;                    % ���͏��

desired = [ 0;0;0;0;0;0; ];

%% C/GMRES�̃R���g���[���p�v�Z
buff = c2d( ss( tf( [ cgmres.tf ], [ (1/cgmres.alpha), 1 ] ) ), dSamplingPeriod );
cgmres.T_outGain = buff.a;              % �\�����ԍ���������
cgmres.T_inGain = buff.b;               % �\�����ԍ���������
clearvars buff;

% �������͒l�̌v�Z�iNewton�@�j
lmd0 = dPhidx( cgmres.x0, desired, cgmres );
cgmres.u0 = [0;0;1;1;0.01;0.01;]; % Newton�@�̏����l

for cnt = 1:20
    cgmres.u0 = cgmres.u0 - ddHddu( cgmres.x0, cgmres.u0, lmd0, sys, cgmres ) \ dHdu( cgmres.x0, cgmres.u0, lmd0, sys, cgmres );
end

cgmres.len_x = length( cgmres.x0 );     % ��Ԃ̐�
cgmres.len_u = length( cgmres.u0 );     % ����ʂ̐�
cgmres.len_lmd = cgmres.len_x;          % �����ϐ��̐�


%% �V�~�����[�V�����̎��s�Ǝ��Ԍv��
tic;
sim( 'CGMRES' );
toc

%% �`��
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

