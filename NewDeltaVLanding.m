%% Landing Site dv Calculations 
% Written by John Smigo
%% Phobos Collection sites
Rp = 11500;
s5 = [145 75];
s18 = [92 57];
s17 = [92 43];
s14 = [175 -4];
s7 = [175 -65];
s16 = [105 -40];
s10 = [37 -70];
s12 = [345 -25];
s4 = [55 -13];
s2 = [52 5];
s3 = [39 6];
s1 = [52 17];
s19 = [39 40];
s20 = [20 19];
s11 = [356 35];
s9 = [315 70];
s21 = [265 57];
s15 = [245 28];
s13 = [295 -5];
s8 = [205 -68];
s6 = [195 78];
let_p = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u"];
path_p = ["ab","bc","cd","de","ef","fg","gh","hi","ij","jk","kl","lm","mn","no","op","pq","qr","rs","st","tu"];
Phobos_Sites_LongLat = struct('a',s6,'b',s5,'c',s14,'d',s8,'e',s7,'f',s16,'g',s10,'h',s12,'i',s4,'j',s2,'k',s3,'l',s1,'m',s17,'n',s18,'o',s19,'p',s20,'q',s11,'r',s9,'s',s21,'t',s15,'u',s13);
for i = 1:21
    j = let_p(i);
    LongLat = Phobos_Sites_LongLat.(j);
    if i < 21
        k = let_p(i+1);
        LongLat2 = Phobos_Sites_LongLat.(k);
        diff = (LongLat2-LongLat)/2;
        rx = Rp*cosd(LongLat(2)+diff(2))*cosd(LongLat(1)+diff(1));
        ry = Rp*cosd(LongLat(2)+diff(2))*sind(LongLat(1)+diff(1));
        rz = Rp*sind(LongLat(2)+diff(2));
        if i == 7
            rx = Rp*cosd(LongLat(2)+diff(2))*cosd(11);
            ry = Rp*cosd(LongLat(2)+diff(2))*sind(11);
            rz = Rp*sind(LongLat(2)+diff(2));
        end
        if i == 8
            rx = Rp*cosd(LongLat(2)+diff(2))*cosd(20);
            ry = Rp*cosd(LongLat(2)+diff(2))*sind(20);
            rz = Rp*sind(LongLat(2)+diff(2));
        end
        if i == 16
            rx = Rp*cosd(LongLat(2)+diff(2))*cosd(8);
            ry = Rp*cosd(LongLat(2)+diff(2))*sind(8);
            rz = Rp*sind(LongLat(2)+diff(2));
        end
        Phobos_Sites_r.(j) = [rx ry rz];
    end
    rx = (Rp-500)*cosd(LongLat(2))*cosd(LongLat(1));
    ry = (Rp-500)*cosd(LongLat(2))*sind(LongLat(1));
    rz = (Rp-500)*sind(LongLat(2));
    Phobos_Sites_gnd.(j) = [rx ry rz];
end
%% Phobos dv Calc
mp = 10.8*10^15;
G = 6.67*10^-11;
mu_p = G*mp;
DM = 1;
dv_p = struct();
ToF_p = struct();
w = 227.25*10^-6;
for z = 2:21
    x = z-1;
    j = let_p(z);
    i = let_p(x);
    k = path_p(x);
    ri_p_up = Phobos_Sites_gnd.(i);
    rf_p_up = Phobos_Sites_r.(i);
    ri_p_down = Phobos_Sites_r.(i);
    rf_p_down = Phobos_Sites_gnd.(j);
    a = 0;
    lastdv_up = 10^5;
    lastdv_down = 10^5;
    for ToF = 200:25:10^4
        a = a+1;
        [v1,v2,C3F] = Lamberts(ri_p_up,rf_p_up,ToF,DM,mu_p);
        [v3,v4,C3F] = Lamberts(ri_p_down,rf_p_down,ToF,DM,mu_p);
        % 2 Norm
        dv1 = norm(v1);
        dv2 = norm(v3-v2);
        dv3 = norm(v4);
        dv_2norm = dv1+dv2+dv3;
        manuev_up = dv1+dv2;
        manuev_down = dv2+dv3;
        % 1 Norm
        dv1norm1 = norm(v1,1);
        dv1norm2 = norm(v3-v2,1);
        dv1norm3 = norm(v4,1);
        manuev_up_pess = dv1norm1+dv1norm2;
        manuev_down_pess = dv1norm2+dv1norm3;
        dv_1norm = dv1norm1+dv1norm2+dv1norm3;
        dvplot_up(a,x) = manuev_up;
        dvplot_down(a,x) = manuev_down;
        if manuev_up < lastdv_up
            v1p(:,x) = v1';
            v2p(:,x) = v2';
            r1p(:,x) = ri_p_up';
            r2p(:,x) = rf_p_up';
            dv_p_up_opt.(k) = manuev_up;
            ToF_p_up.(k) = ToF;
            dv_p_up_pess.(k) = manuev_up_pess;
            dv2o_store.(k) = dv2;
            dv2p_store.(k) = dv1norm2;
            lastdv_up = manuev_up;
        end
        if manuev_down < lastdv_down
            v3p(:,x) = v3';
            v4p(:,x) = v4';
            r3p(:,x) = ri_p_down';
            r4p(:,x) = rf_p_down';
            dv_p_down_opt.(k) = manuev_down;
            ToF_p_down.(k) = ToF;
            dv_p_down_pess.(k) = manuev_down_pess;
            lastdv_down = manuev_down;
        end
        
    end
    dv_p_opt.(k) = dv_p_down_opt.(k)+dv_p_up_opt.(k)-dv2o_store.(k);
    dv_p_pess.(k) = dv_p_down_pess.(k)+dv_p_up_pess.(k)-dv2p_store.(k); 
    ToF_p.(k) = ToF_p_down.(k)+ToF_p_up.(k);
end
figure(1)
hold on
plot([200:25:10^4],dvplot_up(:,1),'r')
plot([200:25:10^4],dvplot_down(:,1),'g')
plot([200:25:10^4],dvplot_down(:,1)+dvplot_up(:,1),'b')
xlabel('ToF (s)')
ylabel('dv (m/s)')
title('dv vs ToF: Phobos Manuever 1')
legend('Up Manuever','Down Manuever','Total dV')
hold off
tot_dv_po = 0;
tot_dv_pp = 0;
tot_ToF_p = 0;
for i = 1:20
    j = path_p(i);
    tot_dv_po = dv_p_down_opt.(j)+dv_p_up_opt.(j)+tot_dv_po-dv2o_store.(j);
    tot_ToF_p = ToF_p_up.(j)+ToF_p_down.(j)+tot_ToF_p;
    tot_dv_pp = dv_p_down_pess.(j)+dv_p_up_pess.(j)+tot_dv_pp-dv2p_store.(k);
end
disp(tot_dv_po)
disp(tot_dv_pp)
disp(tot_ToF_p)
for i = 1:20
    figure(2)
    k = path_p(i);
    [x,y,z] = sphere();
    rad = 11000;
    surf(rad*x,rad*y,rad*z,'FaceColor','none','LineStyle',':')
    hold on
    tup = [0 ToF_p_up.(k)];
    z1 = [r1p(1,i) v1p(1,i) r1p(2,i) v1p(2,i) r1p(3,i) v1p(3,i)];
    tolerance = 10^-13;
    options = odeset('RelTol',tolerance,'AbsTol',tolerance);
    [t1,RV1] = ode45(@propagate_2BP,tup,z1,options,mu_p);
    tdown = [0 ToF_p_down.(k)];
    z2 = [r3p(1,i) v3p(1,i) r3p(2,i) v3p(2,i) r3p(3,i) v3p(3,i)];
    tolerance = 10^-13;
    options = odeset('RelTol',tolerance,'AbsTol',tolerance);
    [t2,RV2] = ode45(@propagate_2BP,tdown,z2,options,mu_p);
    plot3(RV1(:,1),RV1(:,3),RV1(:,5),'r',RV2(:,1),RV2(:,3),RV2(:,5),'r',r1p(1,i),r1p(2,i),r1p(3,i),'b',r2p(1,i),r2p(2,i),r2p(3,i),'gx',r4p(1,i),r4p(2,i),r4p(3,i),'b','LineWidth',2)
    zlabel('z')
end
for i = 1:20
    text(r1p(1,i),r1p(2,i),r1p(3,i),num2str(i),'Color','b','FontWeight','bold','FontSize',12)
end
p(1) = plot3(NaN,NaN,NaN,'r');
p(2) = plot3(NaN,NaN,NaN,'b');
p(3) = plot3(NaN,NaN,NaN,'gx');
text(r4p(1,20),r4p(2,20),r4p(3,20),num2str(21),'Color','b','FontWeight','bold','FontSize',12)
axis equal
title('Phobos Landing Site Trajectories')
xlabel('Inertial x-axis (m)')
ylabel('Inertial y-axis (m)')
zlabel('Inertial z-axis (m)')
legend(p,'Trajectory','# - Landing Site','Midpoint Manuever')
legend('Location','northeast')
%% Deimos Collection sites
Rd = 6700;
s4 = [190 -45];
s14 = [190 -10];
s6 = [207 20];
s2 = [202 80];
s9 = [270 50];
s13 = [260 0];
s15 = [325 -3];
s16 = [355 -55];
s11 = [22 -20];
s17 = [22 10];
s12 = [355 35];
s5 = [35 45];
s8 = [67 72];
s1 = [125 72];
s7 = [125 23];
s10 = [80 -5];
s3 = [115 -45];

let_d = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q"];
path_d = ["ab","bc","cd","de","ef","fg","gh","hi","ij","jk","kl","lm","mn","no","op","pq"];
Deimos_Sites_LongLat = struct('a',s3,'b',s10,'c',s7,'d',s1,'e',s8,'f',s5,'g',s12,'h',s17,'i',s11,'j',s16,'k',s15,'l',s13,'m',s9,'n',s6,'o',s14,'p',s4);
for i = 1:16
    j = let_d(i);
    LongLat = Deimos_Sites_LongLat.(j);
    if i < 16
        k = let_d(i+1);
        LongLat2 = Deimos_Sites_LongLat.(k);
        diff = (LongLat2-LongLat)/2;
        rx = Rd*cosd(LongLat(2)+diff(2))*cosd(LongLat(1)+diff(1));
        ry = Rd*cosd(LongLat(2)+diff(2))*sind(LongLat(1)+diff(1));
        rz = Rd*sind(LongLat(2)+diff(2));
        if i == 6
            rx = Rd*cosd(LongLat(2)+diff(2))*cosd(15);
            ry = Rd*cosd(LongLat(2)+diff(2))*sind(15);
            rz = Rd*sind(LongLat(2)+diff(2));
        end
        if i == 7
            rx = Rd*cosd(LongLat(2)+diff(2))*cosd(9);
            ry = Rd*cosd(LongLat(2)+diff(2))*sind(9);
            rz = Rd*sind(LongLat(2)+diff(2));
        end
        if i == 9
            rx = Rd*cosd(LongLat(2)+diff(2))*cosd(9);
            ry = Rd*cosd(LongLat(2)+diff(2))*sind(9);
            rz = Rd*sind(LongLat(2)+diff(2));
        end
        Deimos_Sites_r.(j) = [rx ry rz];
    end
    rx = (Rd-500)*cosd(LongLat(2))*cosd(LongLat(1));
    ry = (Rd-500)*cosd(LongLat(2))*sind(LongLat(1));
    rz = (Rd-500)*sind(LongLat(2));
    Deimos_Sites_gnd.(j) = [rx ry rz];
end
%% Deimos dv Calc
md = 1.8*10^15;
G = 6.67*10^-11;
mu_d = G*md;
DM = 1;
dv_do = struct();
ToF_d = struct();
w = 57.6*10^-6;
for z = 2:16
    x = z-1;
    j = let_d(z);
    i = let_d(x);
    k = path_d(x);
    ri_d_up = Deimos_Sites_gnd.(i);
    rf_d_up = Deimos_Sites_r.(i);
    ri_d_down = Deimos_Sites_r.(i);
    rf_d_down = Deimos_Sites_gnd.(j);
    b = 0;
    lastdv_up = 10^5;
    lastdv_down = 10^5;
    for ToF = 2*10^2:25:10^4
        b = b+1;
        [v1,v2,C3F] = Lamberts(ri_d_up,rf_d_up,ToF,DM,mu_d);
        [v3,v4,C3F] = Lamberts(ri_d_down,rf_d_down,ToF,DM,mu_d);
        % 2 Norm
        dv1 = norm(v1);
        dv2 = norm(v3-v2);
        dv3 = norm(v4);
        dv_2norm = dv1+dv2+dv3;
        manuev_up = dv1+dv2;
        manuev_down = dv2+dv3;
        % 1 Norm
        dv1norm1 = norm(v1,1);
        dv1norm2 = norm(v3-v2,1);
        dv1norm3 = norm(v4,1);
        manuev_up_pess = dv1norm1+dv1norm2;
        manuev_down_pess = dv1norm2+dv1norm3;
        dv_1norm = dv1norm1+dv1norm2+dv1norm3;
        dvplotd_up(b,x) = manuev_up;
        dvplotd_down(b,x) = manuev_down;
        if manuev_up < lastdv_up
            s = b;
            v1d(:,x) = v1';
            v2d(:,x) = v2';
            r1d(:,x) = ri_d_up';
            r2d(:,x) = rf_d_up';
            dv_d_up_opt.(k) = manuev_up;
            ToF_d_up.(k) = ToF;
            dv_d_up_pess.(k) = manuev_up_pess;
            dv2o_store.(k) = dv2;
            dv2p_store.(k) = dv1norm2;
            lastdv_up = manuev_up;
        end
        if manuev_down < lastdv_down
            q = b;
            v3d(:,x) = v3';
            v4d(:,x) = v4';
            r3d(:,x) = ri_d_down';
            r4d(:,x) = rf_d_down';
            dv_d_down_opt.(k) = manuev_down;
            ToF_d_down.(k) = ToF;
            dv_d_down_pess.(k) = manuev_down_pess;
            lastdv_down = manuev_down;
        end
    end
    dv_d_opt.(k) = dv_d_down_opt.(k)+dv_d_up_opt.(k)-dv2o_store.(k);
    dv_d_pess.(k) = dv_d_down_pess.(k)+dv_d_up_pess.(k)-dv2p_store.(k);  
    ToF_d.(k) = ToF_d_down.(k)+ToF_d_up.(k);
end
figure(3)
hold on
plot([2*10^2:25:10^4],dvplotd_up(:,1),'r')
plot([2*10^2:25:10^4],dvplotd_down(:,1),'g')
plot([200:25:10^4],dvplotd_down(:,1)+dvplotd_up(:,1),'b')
xlabel('ToF (s)')
ylabel('dv (m/s)')
title('dv vs ToF: Deimos Manuever 1')
legend('Up Manuever','Down Manuever','Total dV')
hold off
tot_dv_do = 0;
tot_dv_dp = 0;
tot_ToF_d = 0;
for i = 1:15
    j = path_d(i);
    tot_dv_do = dv_d_down_opt.(j)+dv_d_up_opt.(j)+tot_dv_do-dv2o_store.(j);
    tot_ToF_d = ToF_d_up.(j)+ToF_d_down.(j)+tot_ToF_d;
    tot_dv_dp = dv_d_down_pess.(j)+dv_d_up_pess.(j)+tot_dv_dp;
end
disp(tot_dv_do)
disp(tot_dv_dp)
disp(tot_ToF_d)
for i = 1:15
    figure(4)
    k = path_d(i);
    [x,y,z] = sphere();
    rad = 6200;
    surf(rad*x,rad*y,rad*z,'FaceColor','none','LineStyle',':')
    hold on
    tup = [0 ToF_d_up.(k)];
    z1 = [r1d(1,i) v1d(1,i) r1d(2,i) v1d(2,i) r1d(3,i) v1d(3,i)];
    tolerance = 10^-13;
    options = odeset('RelTol',tolerance,'AbsTol',tolerance);
    [t1,RV1] = ode45(@propagate_2BP,tup,z1,options,mu_d);
    tdown = [0 ToF_d_down.(k)];
    z2 = [r3d(1,i) v3d(1,i) r3d(2,i) v3d(2,i) r3d(3,i) v3d(3,i)];
    tolerance = 10^-13;
    options = odeset('RelTol',tolerance,'AbsTol',tolerance);
    [t2,RV2] = ode45(@propagate_2BP,tdown,z2,options,mu_d);
    plot3(RV1(:,1),RV1(:,3),RV1(:,5),'r',RV2(:,1),RV2(:,3),RV2(:,5),'r',r1d(1,i),r1d(2,i),r1d(3,i),'b',r2d(1,i),r2d(2,i),r2d(3,i),'gx',r4d(1,i),r4d(2,i),r4d(3,i),'b','LineWidth',2)
    zlabel('z')
end
for i = 1:15
    text(r1d(1,i),r1d(2,i),r1d(3,i),num2str(i),'Color','b','FontWeight','bold','FontSize',12)
end
p(1) = plot3(NaN,NaN,NaN,'r');
p(2) = plot3(NaN,NaN,NaN,'b');
p(3) = plot3(NaN,NaN,NaN,'gx');
text(r4d(1,15),r4d(2,15),r4d(3,15),num2str(16),'Color','b','FontWeight','bold','FontSize',12)
axis equal
title('Deimos Landing Site Trajectories')
xlabel('Inertial x-axis (m)')
ylabel('Inertial y-axis (m)')
zlabel('Inertial z-axis (m)')
legend(p,'Trajectory','# - Landing Site','Midpoint Manuever')

function out = propagate_2BP(t,z,mu)
rmag = sqrt(z(1)^2+z(3)^2+z(5)^2);
out(1,1) = z(2);
out(2,1) = -mu/rmag^3*z(1);
out(3,1) = z(4);
out(4,1) = -mu/rmag^3*z(3);
out(5,1) = z(6);
out(6,1) = -mu/rmag^3*z(5);
end
