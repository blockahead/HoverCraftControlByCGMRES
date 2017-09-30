%% ���ݎ�������T�b�����܂ł̐����ϐ��̗\���iEuler�ߎ��j
function lmd = Backward( x, u, d, T, sys, cgmres )
    dt = T / cgmres.dv;
    
    lmd = zeros( cgmres.len_lmd * cgmres.dv, 1 );
    lmd((1:cgmres.len_lmd)+cgmres.len_lmd*(cgmres.dv-1)) = dPhidx( x((1:cgmres.len_x)+cgmres.len_x*(cgmres.dv-1)), d((1:cgmres.len_x)+cgmres.len_x*(cgmres.dv-1)), cgmres );
    
    for cnt = cgmres.dv-1:-1:1
        lmd((1:cgmres.len_lmd)+cgmres.len_lmd*(cnt-1)) = lmd((1:cgmres.len_lmd)+cgmres.len_lmd*(cnt)) ...
                                                            + dHdx( x((1:cgmres.len_x)+cgmres.len_x*(cnt)), u, lmd((1:cgmres.len_lmd)+cgmres.len_lmd*(cnt)), d((1:cgmres.len_x)+cgmres.len_x*(cnt)), sys, cgmres ) * dt;
    end
end