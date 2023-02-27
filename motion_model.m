function [x_post_upd, v_post_upd, R_post_upd, b_w_post_upd, b_a_post_upd, P_post_upd]...
    = motion_model(x_post, v_post, R_post, b_w_post, b_a_post, P_post, w_i, a_i, b_w_dot, ...
                    b_a_dot, Rbi, v_base, Sigma_v, g, T)

w_tilde = Rbi*(w_i - b_w_post);
a_tilde = Rbi*(a_i - b_a_post);

x_pri = x_post + T*v_post;
v_pri = v_post + T*(cross(-w_tilde,v_post) + (R_post)\g + a_tilde);
R_pri = R_post*expm(skew(w_tilde)*T);
b_w_pri = b_w_post + ones(3,1)*b_w_dot*T;
b_a_pri = b_a_post + ones(3,1)*b_a_dot*T;

Gc = [zeros(3,3) R_post         -R_post*skew(v_post);...
      zeros(3,3) -skew(w_tilde) skew(R_post'*g);...
      zeros(3,3) zeros(3,3)     -skew(w_tilde)];
      
Gk = eye(9) + Gc*T;

Vc = [zeros(3,6)
      skew(v_post) eye(3);...
      eye(3)       zeros(3,3)];
  
% This is for MEMS
gyro_var = 0.000025*ones(3,1);
accel_var = 25*(10^-12)*ones(3,1);
Qc = blkdiag(diag(gyro_var), diag(accel_var));

Qk = Vc*Qc*Vc'*T;

P_pri = Gk*P_post*Gk' + Qk;

C = [zeros(3,3) eye(3) zeros(3,3)];

K_gain = (P_pri*C')/(C*P_pri*C' + Sigma_v);

state_post_upd = [x_pri; v_pri; [0;0;0]] + K_gain*(v_base' - v_pri);

x_post_upd = state_post_upd(1:3,1);
v_post_upd = state_post_upd(4:6,1);
R_post_upd = R_pri*expm(skew(state_post_upd(7:9,1)));

P_post_upd = (eye(9) - K_gain*C)*P_pri;

b_w_post_upd = b_w_pri; b_a_post_upd = b_a_pri;

end

