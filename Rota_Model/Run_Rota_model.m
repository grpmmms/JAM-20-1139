%  The Rotational Model
%  -------------------------------------------------------------------------
%  Note:
%      (1) all parameters are dimensionless, see the paper for details
%      (2) the fixed end of the beam defaults to the origin of the coordinate system
%  -------------------------------------------------------------------------
clc
clear
tic

%  m: the scale factor between the object radius and the fixed end of the beam ------{self-definning term}
%  miu ;the interfacial friction coefficient ------{self-definning term}
%  data_p, data_MP: arrays used to store parameters for output
m = 0.2;
miu = 0.1;
data_P = [];
data_MP = [];

%  i: the angle at the point of contact, radian measure ------{self-definning term}
for i = 0.01:0.01:1.1
    theta_m = -atan(miu);
    theta_obj = i;
    
    %  // Call the subfunction Rota_model.m //
    S_p = Rota_model(m,theta_obj,theta_m);
    data_P = [data_P;S_p];
    data_MP = [data_MP;S_p(14),S_p(1),S_p(7),S_p(8),S_p(19),S_p(21)];
end

%  -------------------------------------------------------------------------
%  // Tmage of the axial and tangent force at the beam end along the moving direction //
%  -------------------------------------------------------------------------
theta_x = data_P(:,14);
P_axi_x = data_P(:,19);
P_tri_x = data_P(:,21);
figure
plot(theta_x,P_axi_x)
hold on;
plot(theta_x,P_tri_x)

%  -------------------------------------------------------------------------
%  // Tmage of the contact angle with respect to the angle of the fixed end //
%  -------------------------------------------------------------------------
theta_end = data_P(:,1);
theta_apex = data_P(:,2);
figure
plot(-theta_apex,theta_end,'-ob')

%  -------------------------------------------------------------------------
%  // Tmage of the contact position //
%  -------------------------------------------------------------------------
s_obj = data_P(:,6);
x_obj = data_P(:,7);
y_obj = data_P(:,8);
figure
plot(theta_end, s_obj,'-r')
hold on;
plot(theta_end, x_obj,'-b')
plot(theta_end, y_obj,'or')

%  -------------------------------------------------------------------------
%  // Another image of the contact position to compare with the former //
%  -------------------------------------------------------------------------
x_obj_1 = data_P(:,9);
y_obj_1 = data_P(:,10);
hold on;
plot(theta_end, x_obj_1,'ob')
plot(theta_end, y_obj_1,'-r')

%  -------------------------------------------------------------------------
%  // Tmage of the X and Y component force at the beam end along the moving direction //
%  -------------------------------------------------------------------------
F_obj = data_P(:,11);
Fx_obj = data_P(:,12);
Fy_obj = data_P(:,13);
figure
plot(theta_end, Fy_obj,'-b')
hold on;
plot(theta_end, Fx_obj,'ob')

%  -------------------------------------------------------------------------
%  // Tmage of the bending moment at the beam end along the moving direction //
%  -------------------------------------------------------------------------
theta_M_apex = data_P(:,14);
M_apex = data_P(:,15);
figure
plot(theta_M_apex,M_apex,'ob')

%  -------------------------------------------------------------------------
%  // Tmage of the strain energy of beam along the moving direction //
%  -------------------------------------------------------------------------
U_apex = data_P(:,16);
figure
plot(theta_M_apex,U_apex,'or')

%  -------------------------------------------------------------------------
%  // Tmage of the derivative of the strain energy along the moving direction //
%  -------------------------------------------------------------------------
%  Output the derivative of the strain energy 
z = gradient(U_apex)./gradient(theta_M_apex);
figure
plot(theta_M_apex,z,'or')
t = gradient(z)./gradient(theta_M_apex);
figure
plot(theta_M_apex,t,'ob')
z_t = [];
[row,col] = size(z);
for i = 1:1:row
    z_1 = z(i,1);
    z_2 = t(i,1);
    theta_aa = theta_M_apex(i,1);
    theta_u = U_apex(i,1);
    z_t = [z_t;z_1 z_2];
end
R_es = z_t;
fid = fopen('data__U_dev.txt','w');
fprintf(fid,'%g\n');
[row,col] = size(R_es);
for i = 1:1:row
    for j = 1:1:col
        if (j == col)
            fprintf(fid,'%g\n',R_es(i,j));
        else
            fprintf(fid,'%g\t',R_es(i,j));
        end
    end
end
fclose(fid);

%  -------------------------------------------------------------------------
%  // If necessary, output the parameters you need. //
%  -------------------------------------------------------------------------
R_es = data_MP;
fid = fopen('Total-Parameters_Rota.txt','w');
fprintf(fid,'%g\n');
[row,col] = size(R_es);
for i = 1:1:row
    for j = 1:1:col
        if (j == col)
            fprintf(fid,'%g\n',R_es(i,j));
        else
            fprintf(fid,'%g\t',R_es(i,j));
        end
    end
end
fclose(fid);
toc