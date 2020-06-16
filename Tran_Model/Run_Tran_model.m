%  The Translation Model
%  -------------------------------------------------------------------------
%  Note:
%      (1) all parameters are dimensionless, see the paper for details
%      (2) the fixed end of the beam defaults to the origin of the coordinate system
%  -------------------------------------------------------------------------
clc
clear
tic

%  F_F: an array used to store forces 
%  S_F: an array used to store location coordinate position 
%  ...Arrays used to store parameters for output
F_F = [];
S_F = [];
F_S = [];
F_totel = [];
F_to_1 = [];
F_to_2 = [];

%  m: the scale factor between the object radius and the fixed end of the beam -------{self-definning term}
%  miu ;the interfacial friction coefficient ------{self-definning term}
m = 0.1;
miu = 0;

%  i: the angle at the point of contact, radian measure ------{self-definning term}
for i = 0.01:0.02:1.5
    fhib = i;
    fhia = -atan(miu);
    
    %  Integral term in arc length and coordinate formula (Eq.7) int the paper
    fun1 = @(s)cos(s)./sqrt((cos(fhib+fhia)*sin(fhib)-sin(fhib+fhia)*cos(fhib))+(sin(fhib+fhia).*cos(s)-cos(fhib+fhia).*sin(s)));
    fun2 = @(s)sin(s)./sqrt((cos(fhib+fhia)*sin(fhib)-sin(fhib+fhia)*cos(fhib))+(sin(fhib+fhia).*cos(s)-cos(fhib+fhia).*sin(s)));
    fun3 = @(s)1./sqrt((cos(fhib+fhia)*sin(fhib)-sin(fhib+fhia)*cos(fhib))+(sin(fhib+fhia).*cos(s)-cos(fhib+fhia).*sin(s)));
    F_mid = integral(fun1,0,fhib,'ArrayValued',true);
    F_mid_1=integral(fun2,0,fhib,'ArrayValued',true);
    F_mid_2=integral(fun3,0,fhib,'ArrayValued',true);
    kt = 1./(1-m*sin(fhib));
       
    %  F is equal to 2*sqrt(alpha), where alpha is a dimensionless parameter of force.
    F = kt*F_mid;
    
    %  F_p: the resultant force
    %  F_i: the component of the resultant force in the Y direction, its value can be positive or negative
    %  F_j: the component of the resultant force in the X direction, its value can be only positive
    
    F_p = 0.5.*F.^2; 
    F_i = F_p.*cos(fhib+fhia);
    F_j = F_p.*sin(fhib+fhia);
    F_F = [F_F;i F_i F_j F_p];
    
    %  s_obj: the arc length at the point of contact
    %  y_obj: the y-coordinate of the contact point
    %  x_obj: the x-coordinate of the contact point
    s_obj = F_mid_2/F;
    y_obj = F_mid_1/F;
    x_obj = 1-m*sin(fhib);
    S_F = [S_F;i s_obj x_obj y_obj];
   
    %  Y_true:the displacement of object center in the y direction
    %  Y_tran: The moving distance of the fixed end of the beam in the y direction
    Y_true = y_obj-m*cos(fhib);  
    Y_tran = y_obj-m*cos(fhib)+m;
    
    %  -------------------------------------------------------------------------
    %  // The process to solve the strain energy of the beam //
    %  -------------------------------------------------------------------------
    s_alpha = 0.5.*F;
    theta_s = 0;
    theta_obj = fhib;
    theta_int = theta_obj+fhia;
    
    %  M_i ; an array used to store bending moments temporarily
    M_i = [];
    d_theta = (theta_obj-theta_s)/200;
    for theta_ti = theta_s:d_theta:theta_obj-d_theta
        funa = @(s)cos(s)./((cos(theta_int).*sin(theta_obj)-sin(theta_int).*cos(theta_obj))+(sin(theta_int).*cos(s)-cos(theta_int).*sin(s))).^0.5;
        funb = @(s)sin(s)./((cos(theta_int).*sin(theta_obj)-sin(theta_int).*cos(theta_obj))+(sin(theta_int).*cos(s)-cos(theta_int).*sin(s))).^0.5;
        aa = theta_s;
        bb = theta_ti;
        inta = integral(funa,aa,bb,'ArrayValued',true);
        intb = integral(funb,aa,bb,'ArrayValued',true);
        
        M_ti = ((1-m*sin(theta_obj)-1./2./s_alpha.*inta).*cos(theta_int))+((Y_true+m.*cos(theta_obj))-1./2./s_alpha.*intb).*sin(theta_int);
        M_t_n = M_ti.*(2.*(s_alpha).^2);  
        int_s = 1./((cos(theta_int).*sin(theta_obj)-sin(theta_int).*cos(theta_obj))+(sin(theta_int).*cos(theta_ti)-cos(theta_int).*sin(theta_ti))).^0.5;
        M_tj = 1/4./s_alpha.*int_s.*M_t_n.^2;
        M_i = [M_i;M_tj];
    end
    
    %  U_tj : the strain energy of the beam
    theta_tj = theta_s:d_theta:theta_obj-d_theta;
    theta_tj = theta_tj';
    U_tj = trapz(theta_tj,M_i);
    
    %  -------------------------------------------------------------------------
    %  // Estimate the beam stability //
    %  -------------------------------------------------------------------------
    Fx_obj = F_j;
    Fy_obj = F_i;
    S_obj = s_obj;
    x_0 = 0;
    y_0 = 0;
    gama = atan((-y_obj-y_0)./(x_obj-x_0));
    
    %  P_axi: the axial force along the beam
    %  P_txi: the tangent force along the beam
    %  P_c: the bucking critical force
    P_axi = Fx_obj.*cos(gama)+Fy_obj.*sin(gama);
    P_axi = P_axi.*S_obj.^2;
    P_tri = -Fx_obj.*sin(gama)+Fy_obj.*cos(gama);
    P_tri = P_tri.*S_obj.^2;
    P_c = pi.*pi;
    
    F_S = [F_S;i Fx_obj Fy_obj F_p i s_obj x_obj y_obj];
    F_totel = [F_totel;Y_tran Fy_obj Fx_obj];
    F_to_1 = [F_to_1;Y_tran U_tj s_obj fhib Fy_obj];
    F_to_2 = [F_to_2;i Y_tran s_obj P_axi P_tri s_alpha];

end

%  -------------------------------------------------------------------------
%  // Tmage of the axial and tangent force at the beam end along the moving direction //
%  -------------------------------------------------------------------------
Y_tran_x = F_to_2(:,2);
P_axi_x = F_to_2(:,3);
P_tri_x = F_to_2(:,4);
figure
plot(Y_tran_x,P_axi_x,'ob')
plot(Y_tran_x,P_tri_x,'-b')

%  -------------------------------------------------------------------------
%  // Tmage of the X and Y component force at the beam end along the moving direction //
%  -------------------------------------------------------------------------
theta = F_F(:,1);
Ft = F_F(:,2);
Fn = F_F(:,3);
Fp = F_F(:,4);
figure
plot(theta, Ft,'-b')
hold on;
plot(theta, Fn,'ob')

%  -------------------------------------------------------------------------
%  // Tmage of the contact position //
%  -------------------------------------------------------------------------
St = S_F(:,2);
Sx = S_F(:,3);
Sy = S_F(:,4);
figure
plot(theta, St,'-r')
hold on;
plot(theta, Sx,'-r')
hold on;
plot(theta, Sy,'or')

%  -------------------------------------------------------------------------
%  // Tmage of the X and Y component force at the beam end along the relative distance //
%  -------------------------------------------------------------------------
Sy_true = F_totel(:,1);
Fy_true = F_totel(:,2);
Fx_true = F_totel(:,3);
figure
plot(Sy_true, Fy_true,'or')
hold on;
plot(Sy_true, Fx_true,'ob')

%  -------------------------------------------------------------------------
%  // Tmage of the strain energy of beam along the moving direction //
%  -------------------------------------------------------------------------
Sy_tran = F_to_1(:,1);
U_tran = F_to_1(:,2);
figure
plot(Sy_tran, U_tran,'or')

%  -------------------------------------------------------------------------
%  // If necessary, the load ratio during contact can be calaulated. //
%  -------------------------------------------------------------------------
[row,colum] = size(F_totel);
tt = [0];
ft = [0];
fn = [0];
sumU = [];
for w = 1:row
    Sy_y = F_totel(w,1);
    Fy_y = F_totel(w,2);
    Fx_x = F_totel(w,3);
    tt = [tt;Sy_y];
    ft = [ft;Fy_y];
    fn = [fn;Fx_x];
    tt = tt';
    ft = ft';
    fn = fn';
    sumy = trapz(tt,ft);
    sumx = trapz(tt,fn);
    sunf = sumy/sumx;
    sumU = [sumU;sunf];
    tt = tt';
    ft = ft';
    fn = fn';
end

a=sumU;
fid = fopen('f-load-ratio.txt','wt');
[row,col] = size(a);
for i = 1:1:row
    for j = 1:1:col
        if (j == col)
            fprintf(fid,'%g\n',a(i,j));
        else
            fprintf(fid,'%g\t',a(i,j));
        end
    end
end
fclose(fid);

%  -------------------------------------------------------------------------
%  // If necessary, output the parameters you need. //
%  -------------------------------------------------------------------------
R_es = F_to_2;
fid = fopen('Total-Parameters_Tran.txt','wt');
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