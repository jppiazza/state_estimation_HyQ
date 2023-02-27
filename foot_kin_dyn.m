function [v_base_foot, f_foot] = foot_kin_dyn(foot, sensors, index, Rbi, b_w, l0, l1, l2, l3, dlr, dfh)

switch foot
    case "LF"
        q0 = sensors.LF_HAA_q(index); q0_dot = sensors.LF_HAA_qd_f(index);
        q1 = sensors.LF_HFE_q(index); q1_dot = sensors.LF_HFE_qd_f(index);
        q2 = sensors.LF_KFE_q(index); q2_dot = sensors.LF_KFE_qd_f(index);
        tau_0 = sensors.LF_HAA_forceMeasured(index);
        tau_1 = sensors.LF_HFE_forceMeasured(index)*l1; 
        tau_2 = sensors.LF_KFE_forceMeasured(index)*l2;
    case "RF"
        q0 = sensors.RF_HAA_q(index); q0_dot = sensors.RF_HAA_qd_f(index);
        q1 = sensors.RF_HFE_q(index); q1_dot = sensors.RF_HFE_qd_f(index);
        q2 = sensors.RF_KFE_q(index); q2_dot = sensors.RF_KFE_qd_f(index);
        tau_0 = sensors.RF_HAA_forceMeasured(index);
        tau_1 = sensors.RF_HFE_forceMeasured(index)*l1;
        tau_2 = sensors.RF_KFE_forceMeasured(index)*l2;
    case "LH"
        q0 = sensors.LH_HAA_q(index); q0_dot = sensors.LH_HAA_qd_f(index);
        q1 = sensors.LH_HFE_q(index); q1_dot = sensors.LH_HFE_qd_f(index);
        q2 = sensors.LH_KFE_q(index); q2_dot = sensors.LH_KFE_qd_f(index);
        tau_0 = sensors.LH_HAA_forceMeasured(index);
        tau_1 = sensors.LH_HFE_forceMeasured(index)*l1;
        tau_2 = sensors.LH_KFE_forceMeasured(index)*l2;
    case "RH"
        q0 = sensors.RH_HAA_q(index); q0_dot = sensors.RH_HAA_qd_f(index);
        q1 = sensors.RH_HFE_q(index); q1_dot = sensors.RH_HFE_qd_f(index);
        q2 = sensors.RH_KFE_q(index); q2_dot = sensors.RH_KFE_qd_f(index);
        tau_0 = sensors.RH_HAA_forceMeasured(index);
        tau_1 = sensors.RH_HFE_forceMeasured(index)*l1;
        tau_2 = sensors.RH_KFE_forceMeasured(index)*l2;
end

x_foot_leg = [-l1*sin(q1)-l2*sin(q1+q2);tan(q0)*(-l0-l1*cos(q1)-l2*cos(q1+q2)-l3);-l0-l1*cos(q1)-l2*cos(q1+q2)-l3];

J = [0 -l1*cos(q1)-l2*cos(q1 + q2) -l2*cos(q1 + q2);
    (sec(q0))^2*(-l0-l1*cos(q1)-l2*cos(q1+q2)-l3) tan(q0)*(l1*sin(q1)+l2*sin(q1 + q2)) tan(q0)*(l2*sin(q1 + q2));
    0 l1*sin(q1)+l2*sin(q1 + q2) l2*sin(q1 + q2)];

v_foot = J*[q0_dot; q1_dot; q2_dot];  

% Convert the foot position to the base frame
switch foot
    case "LF"
        x_foot_base = x_foot_leg + [dfh/2; dlr/2; 0];       
    case "RF"
        x_foot_leg(2) = -x_foot_leg(2);
        x_foot_base = x_foot_leg + [dfh/2; -dlr/2; 0];
        v_foot(2) = -v_foot(2);
    case "LH"
        x_foot_base = x_foot_leg + [-dfh/2; dlr/2; 0];
    case "RH"
        x_foot_leg(2) = -x_foot_leg(2);
        x_foot_base = x_foot_leg + [-dfh/2; -dlr/2; 0]; 
        v_foot(2) = -v_foot(2);
end    

w_imu_x = sensors.B_Ad_A_r(index);
w_imu_y = sensors.B_Ad_B_r(index);
w_imu_z = sensors.B_Ad_G_r(index);

w_base_base = Rbi*([w_imu_x; w_imu_y; w_imu_z] - b_w);

v_base_foot = -v_foot - cross(w_base_base, x_foot_base);

Tau = [tau_0; tau_1; tau_2];

% Ground reaction force
f_foot = -J'\Tau;

end

