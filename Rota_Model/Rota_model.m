function S_p = Rota_model(m,theta_obj,theta_m)
theta_int = theta_obj+theta_m;
[x,fval,exitflag] = fsolve(@S_theta,[-1,1]',optimset('Display','off'));
    
    %  // The subfunction to solve the equation set //
    %  one unknown term {x(1)} is the angle at the fixed end,
    %  and the other {x(2)} is a dimensionless parameter of force.
    function y = S_theta(x)
        theta_s = x(1); 
        s_alpha = x(2);
        
        int1_a = sqrt((cos(theta_int).*sin(theta_obj)-sin(theta_int).*cos(theta_obj))+(sin(theta_int).*cos(theta_s)-cos(theta_int).*sin(theta_s)));
        int1_b = (1-m.*sin(theta_obj)).*cos(theta_int)-m.*(1-cos(theta_obj)).*sin(theta_int);
        int1 = int1_a./int1_b;
        y(1) = int1-s_alpha; 
        
        fun = @(s)cos(s)./((cos(theta_int).*sin(theta_obj)-sin(theta_int).*cos(theta_obj))+(sin(theta_int).*cos(s)-cos(theta_int).*sin(s))).^0.5;
        a = theta_s;
        b = theta_obj;
        int2 = integral(fun,a,b,'ArrayValued',true);
        int2_a = 1./2./(1-m.*sin(theta_obj));
        y(2) = int2_a.*int2-s_alpha; 
    end

%  M_i ; an array used to store bending moments temporarily
M_i = [];
d_theta = (theta_obj-theta_s)/100;
for theta_ti = theta_s:d_theta:theta_obj
    funa = @(s)cos(s)./((cos(theta_int).*sin(theta_obj)-sin(theta_int).*cos(theta_obj))+(sin(theta_int).*cos(s)-cos(theta_int).*sin(s))).^0.5;
    funb = @(s)sin(s)./((cos(theta_int).*sin(theta_obj)-sin(theta_int).*cos(theta_obj))+(sin(theta_int).*cos(s)-cos(theta_int).*sin(s))).^0.5;
    aa = theta_s;
    bb = theta_ti;
    inta = integral(funa,aa,bb,'ArrayValued',true);
    intb = integral(funb,aa,bb,'ArrayValued',true);

    M_ti = ((1-m*sin(theta_obj)-1./2./s_alpha.*inta).*cos(theta_int))+((cos(theta_obj)-1).*m-1./2./s_alpha.*intb).*sin(theta_int);
    M_t_n = M_ti.*(2.*(s_alpha).^2);  %normalizin M
    int_s = 1./((cos(theta_int).*sin(theta_obj)-sin(theta_int).*cos(theta_obj))+(sin(theta_int).*cos(theta_ti)-cos(theta_int).*sin(theta_ti))).^0.5;
    M_tj = 1/4./s_alpha.*int_s.*M_t_n.^2;    %normalizin U
    M_i = [M_i;M_tj];
end

%  U_tj : the strain energy of the beam
theta_tj = theta_s:d_theta:theta_obj;
U_tj = trapz(theta_tj,M_i);

%  // Calculate the position of contact //
%  Integral term in arc length and coordinate formula (Eq.7) int the paper
%  s_obj: the arc length at the point of contact
%  x_obj: the x-coordinate of the contact point
%  y_obj: the y-coordinate of the contact point
%  d_obj: the distance between the beam fixed end and the contact point
fun1 = @(s)1./((cos(theta_int).*sin(theta_obj)-sin(theta_int).*cos(theta_obj))+(sin(theta_int).*cos(s)-cos(theta_int).*sin(s))).^0.5;
fun2 = @(s)cos(s)./((cos(theta_int).*sin(theta_obj)-sin(theta_int).*cos(theta_obj))+(sin(theta_int).*cos(s)-cos(theta_int).*sin(s))).^0.5;
fun3 = @(s)sin(s)./((cos(theta_int).*sin(theta_obj)-sin(theta_int).*cos(theta_obj))+(sin(theta_int).*cos(s)-cos(theta_int).*sin(s))).^0.5;
a1 = x(1);
b1 = theta_obj;
int_s = integral(fun1,a1,b1,'ArrayValued',true);
int_x = integral(fun2,a1,b1,'ArrayValued',true);
int_y = integral(fun3,a1,b1,'ArrayValued',true);

para_a = 1./2./x(2);
s_obj = int_s.*para_a;
x_obj = int_x.*para_a;
y_obj = int_y.*para_a;
d_obj = (x_obj.^2+y_obj.^2).^0.5;

%  // Another way to calculate the position of contact //
x_obj_1 = 1-m.*sin(theta_obj);
y_obj_1 = (cos(theta_obj)-1).*m;

%  F_obj : the resultant force
%  Fx_obj : the component of the resultant force in the X direction, its value can be positive or negative
%  Fy_obj : the component of the resultant force in the Y direction, its value can be only positive
F_obj = 0.5.*(2.*x(2)).^2;
Fx_obj = sin(theta_int).*F_obj;
Fy_obj = cos(theta_int).*F_obj;

theta_e = theta_obj;
theta_a = x(1);
theta_c = theta_m;

%  M_apex: the bending moment
int_m = (1-m.*sin(theta_obj)).*cos(theta_int)-m.*(1-cos(theta_obj)).*sin(theta_int);
M_apex = 2*int_m*x(2).^2;

%  -------------------------------------------------------------------------
%  // Estimate the beam stability //
%  -------------------------------------------------------------------------
x_0 = 0;
y_0 = 0;
gama = atan((-y_obj-y_0)./(x_obj-x_0));

%  P_axi: the axial force along the beam
%  P_txi: the tangent force along the beam
%  P_c: the bucking critical force
P_axi = Fx_obj.*cos(gama)+Fy_obj.*sin(gama);
P_axi = P_axi.*s_obj.^2;
P_tri = -Fx_obj.*sin(gama)+Fy_obj.*cos(gama);
P_tri = P_tri.*s_obj.^2;
P_c = pi.*pi;

%  // Output the parameters you need. //
S_p = [theta_e,theta_a,theta_c,x(2),exitflag,s_obj,x_obj,y_obj,x_obj_1,y_obj_1,F_obj,Fx_obj,Fy_obj,-x(1),M_apex,U_tj,d_obj,gama,P_axi,P_c,P_tri,s_alpha];

end
