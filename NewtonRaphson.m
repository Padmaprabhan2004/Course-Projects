clear;
clc;

f_tar=[0,0]';
vec_old=[1,1]';
J=jacobian(vec_old);
vec=vec_old+pinv(J)*(f_tar-function_at_vec(vec_old))
for i=1:1
    vec_old=vec;
    J=jacobian(vec_old);
    vec=vec_old+pinv(J)*(f_tar-function_at_vec(vec_old))
end

function J=jacobian(vec)
    J=[2*vec(1),0;
        0,2*vec(2)];
end

function f=function_at_vec(vec)
    f=[vec(1)^2-9,vec(2)^2-4]';
end