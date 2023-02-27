function [x, v, R, b_w, b_a, P, f_feet, v_base] = HyQ(sensors, x0, v0, R0, b_w_0, b_a_0, P0, IMU)

T = 0.001;
K = size(sensors.time, 1);

g = [0 0 -9.80655]'; % gravity in the world frame
 
% HyQ parameters, all in meters
l0 = 0.08;
l1 = 0.35;
l2 = 0.346;
l3 = 0.02;
dlr = 0.414;
dfh = 0.747;

% Rotation matrices for the IMUs
if IMU == "GX5"
    Rbi = [1 0 0;0 -1 0;0 0 -1];
else
    Rbi = [0 1 0;-1 0 0;0 0 -1];
end

% Bias rates
b_w_dot = 0.0022; % deg/s
b_a_dot = 0.04*10^-3; % g

feet = ["RF","LF","RH","LH"];

f_feet = NaN(K,3,4);
c_feet = zeros(K,4);
v_base_feet = zeros(K,3,4);
v_base = zeros(K,3);

x = NaN(K+1,3); v = NaN(K+1,3); R = NaN(3*(K+1),3);
b_w = NaN(K+1,3); b_a = NaN(K+1,3); P = NaN(9*(K+1),9);

x_hat = x0; v_hat = v0; R_hat = R0; 
b_w_hat = b_w_0; b_a_hat = b_a_0; P_hat = P0;

x(1,:) = x0'; v(1,:) = v0'; R(1:3,1:3) = R0; 
b_w(1,:) = b_w_0'; b_a(1,:) = b_a_0'; P(1:9,1:9) = P0;

var = [0.3 0.3 0.3]; % Values can change
Sigma_0 = diag(var);

for i = 1:K-1
    for j = 1:4
        [v_base_foot, f_foot] = foot_kin_dyn(feet(j), sensors, i, Rbi, b_w_hat, l0, l1, l2, l3, dlr, dfh);
        v_base_feet(i,:,j) = v_base_foot;
        f_feet(i,:,j) = f_foot;
        if f_feet(i,3,j) > 2000
            c_feet(i,j) = 1;
        end
    end
    feet_in_stance = sum(c_feet(i,:));
    if feet_in_stance == 0
        v_base(i,:) = zeros(1,3);
        Sigma_v = Sigma_0;
    else
        v_base(i,:) = [c_feet(i,:)*reshape(v_base_feet(i,1,:),[4,1]);...
                       c_feet(i,:)*reshape(v_base_feet(i,2,:),[4,1]);...
                       c_feet(i,:)*reshape(v_base_feet(i,3,:),[4,1])]*(1/feet_in_stance);
        Var = reshape((v_base_feet(i,:,:) - repmat(v_base(i,:),[1,1,4])),[3,4]);   
        Sigma_m = zeros(3,3);
        for k = 1:4
            Sigma_m = Sigma_m + c_feet(i,k)*(Var(:,k)*Var(:,k)');
        end
        Sigma_v = Sigma_0 + (1/feet_in_stance)*Sigma_m;
    end
   
    if IMU == "GX5"
        w_i = [sensors.B_Ad_A_r(i); sensors.B_Ad_B_r(i); sensors.B_Ad_G_r(i)];
        a_i = [sensors.B_Xdd_r(i); sensors.B_Ydd_r(i); sensors.B_Zdd_r(i)];
    else
        w_i = [sensors.B2_Ad_A_r(i); sensors.B2_Ad_B_r(i); sensors.B2_Ad_G_r(i)];
        a_i = [sensors.B2_Xdd_r(i); sensors.B2_Ydd_r(i); sensors.B2_Zdd_r(i)];
    end
    
    [x_hat_upd, v_hat_upd, R_hat_upd, b_w_hat_upd, b_a_hat_upd, P_hat_upd] = ...
         motion_model(x_hat, v_hat, R_hat, b_w_hat, b_a_hat, P_hat, w_i, a_i, b_w_dot, ...
          b_a_dot, Rbi, v_base(i,:), Sigma_v, g, T);
      
    x(i+1,:) = x_hat_upd'; v(i+1,:) = v_hat_upd'; R(3*i+1:3*i+3,1:3) = R_hat_upd;
    b_w(i+1,:) = b_w_hat_upd'; 
    b_a(i+1,:) = b_a_hat_upd'; 
    P(9*i+1:9*i+9,1:9) = P_hat_upd;
    
    x_hat = x_hat_upd;
    v_hat = v_hat_upd;
    R_hat = R_hat_upd;
    b_w_hat = b_w_hat_upd;
    b_a_hat = b_a_hat_upd;
    P_hat = P_hat_upd;
       
end

end

