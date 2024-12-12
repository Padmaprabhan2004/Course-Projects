%% Initiation
% to type- matrixlog6,adjoint,matrixexp6
T_sd = [1, 0, 0, 0.5;
        0, 1, 0, 0;
        0, 0, 1, 0.4;
        0, 0, 0, 1];
% WAM Dimensions
l1 = 550;
l2 = 300;
l3 = 60;
w1 = 45;

M = [1, 0, 0, 0;
     0, 1, 0, 0;
     0, 0, 1, l1+l2+l3;
     0, 0, 0, 1];

S = [0, 0, 0, 0, 0, 0, 0;
     0, 1, 0, 1, 0, 1, 0;
     1, 0, 1, 0, 1, 0, 1;
     0, l1+l2+l3, 0, l2+l3, 0, l3, 0;
     0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, w1, 0, 0, 0];

% Initial guess
theta = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]';
epsilon_v = 0.0001;
epsilon_w = 0.001;
T_sb = forward_kinematics_space(M, S, theta);
v_s_matrix = MatrixLog6(T_sb \ T_sd);
v_s = se3_to_vec(v_s_matrix);

%% Iterations & plotting
while (norm(v_s(1:3)) > epsilon_w) || (norm(v_s(4:6)) > epsilon_v)
    theta = theta + pinv(jacobian_space(S, theta)) * v_s;
    T_sb = forward_kinematics_space(M, S, theta);
    v_s_matrix = MatrixLog6(T_sb \ T_sd);
    v_s = se3_to_vec(v_s_matrix)
    pause(0.2);
end


  
%% Helper functions

function T=forward_kinematics_space(M,S,theta)
    T=eye(4,4);
    for i=1:length(theta)
        T=T*MatrixExp6(vec_to_se3(S(:,i))*theta(i));
    end
    T=T*M;
end

function se3=vec_to_se3(vec)
    se3=[0,-vec(3),vec(2),vec(4);
        vec(3),0,-vec(1),vec(5);
        -vec(2),vec(1),0,vec(6);
        0,0,0,0];
end

function J_s=jacobian_space(S,theta)
    J_s=S;
    T=eye(4,4);
    for i=2:length(theta)
        T=T*MatrixExp6(vec_to_se3(S(:,i-1))*theta(i-1));
        J_s(:,i)=Adjoint(T)*S(:,i);
    end
end

function v=se3_to_vec(se3)
    v=[se3(3,2);se3(1,3);se3(2,1);se3(1:3,4)];
end


