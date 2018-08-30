%Saimanoj Siddula
%AerE 451 Astrodynamics II
%Midterm Problem 2/3
clear,clc
close all

%times = [11+(20/60),11+(35/60),12]; %hours
%alpha = [-54.1249012,-19.1381645,98.7739537]; %degrees
%delta = [-26.6156514,3.9172775,31.1314513]; %degrees
times = [11+(20/60),11+(30/60),11+(40/60)]; %hours
alpha = [-54.1249012,-33.0588410,-00.4172870]; %degrees
delta = [-22.6156514,-7.2056382,17.4626616]; %degrees

Y= 2007;
M=08;
D=20;

YT = 2000;

%Enter East Longitude angle in degrees
%EL = -110;
EL = 250;

%Finding the Julian Day Given UT---------------------------
J01 = 367*Y;
J02 = (Y+floor((M+9)/12));
J03 = floor(7/4*J02);
J04 = floor((275*M)/9);
J05 = D + 1721013.5;

J0 = J01-J03+J04+J05;

%Julian Date at year 2000
J01T = 367*YT;
J02T = (YT+floor((1+9)/12));
J03T = floor(7/4*J02T);
J04T = floor((275*1)/9);
J05T = 1 + 1721014;

%J0 is Julian Date at UT = 0
JT = J01T-J03T+J04T+J05T;

J01D = 367*Y;
J02D = (Y+floor((1+9)/12));
J03D = floor(7/4*J02D);
J04D = floor((275*1)/9);
J05D = 1 + 1721013.5;

%J0 is the Julian Date at UT =0
J = J01D-J03D+J04D+J05D;
fprintf('The Julian Date at UT is: %.13f\n',J)

%To find Julian Date values after 2000
T0 = (J0-JT)/36525;
DAY = J0-J

%Finding Greenwhich Sidereal Time
G0 = 100.4606184+36000.77004*T0+.000387933*(T0^2)-(2.583*10^-8)*(T0^3);
G01 = G0;

n = 1;
while n==1 %Resetting answer to angle between 360 and 0.
   if G0<0
      G0 = G0+360;
   elseif G0>360
       G0 = G0-360;
   elseif 0<=G0<=360
       n=0;       
   end    
end

GD = G0 + 360.98564724.*times./24
GD1 = GD;

n=1;
while n==1
   if GD<0
       GD = GD+360;
   elseif GD >360
       GD = GD-360;
   elseif 0<=GD<=360
       n=0;       
   end    
end


GR = GD.*pi/180;

%finds local sidereal time
LSTD = GD+EL;
fprintf('The Local Sidereal Day is %.13f\n',LSTD)
LSTR = LSTD.*pi/180;

%Defines Constants for the problem
L = 40;
ae = 6378.1;
ee = .08182;
mu = 3.986e5;
t = (times-(11+(20/60)))*60*60; %time in seconds
tau_1 = t(1)-t(2);
tau_3 = t(3)-t(2);
tau = tau_3 - tau_1;
theta = LSTD;
H = 2.0;

X_val = (ae/(sqrt(1-ee^2*sind(L)^2))+H)*cosd(L);
Z = ((ae*(1-ee^2))/(sqrt(1-ee^2*sind(L)^2))+H)*sind(L);

L_1 = [cosd(alpha(1))*cosd(delta(1));sind(alpha(1))*cosd(delta(1));sind(delta(1))];
L_2 = [cosd(alpha(2))*cosd(delta(2));sind(alpha(2))*cosd(delta(2));sind(delta(2))];
L_3 = [cosd(alpha(3))*cosd(delta(3));sind(alpha(3))*cosd(delta(3));sind(delta(3))];
fprintf('The L(slant direction) values for time 1 is: [%13f,%13f,%13f][km] \n',L_1)
fprintf('The L(slant direction) values for time 1 is: [%13f,%13f,%13f][km] \n',L_2)
fprintf('The L(slant direction) values for time 1 is: [%13f,%13f,%13f][km] \n',L_3)
R_1 = [X_val*cosd(theta(1));X_val*sind(theta(1));Z];
R_2 = [X_val*cosd(theta(2));X_val*sind(theta(2));Z];
R_3 = [X_val*cosd(theta(3));X_val*sind(theta(3));Z];
fprintf('The R values for time 1 is: [%13f;%13f;%13f][km] \n',R_1)
fprintf('The R values for time 2 is: [%13f;%13f;%13f][km] \n',R_2)
fprintf('The R values for time 3 is: [%13f;%13f;%13f][km] \n',R_3)
D_0      = dot(L_1,(cross(L_2,L_3)));                                    %tripple scalar product (rho)
D_11     = dot(R_1, (cross(L_2, L_3)));
D_21     = dot(R_2, (cross(L_2, L_3)));
D_31     = dot(R_3, (cross(L_2, L_3)));
D_12     = dot(R_1,(cross(L_1, L_3)));
D_22     = dot(R_2,(cross(L_1, L_3)));
D_32     = dot(R_3,(cross(L_1, L_3)));
D_13     = dot(R_1,(cross(L_1, L_2)));
D_23     = dot(R_2,(cross(L_1, L_2)));
D_33     = dot(R_3,(cross(L_1, L_2)));

A       = (1/D_0)*(-D_12*(tau_3/tau)+D_22+D_32*(tau_1/tau));
B       = (1/(6*D_0))*(D_12*(tau_3^2-tau^2)*(tau_3/tau)+D_32*(tau^2-tau_1^2)*(tau_1/tau));
E       = dot(R_2, L_2);
fprintf('The A value is %.13f\n',A)
fprintf('The B value is %.13f\n',B)
fprintf('The E value is %.13f\n',E)

a       = -(A^2+2*A*E+norm(R_2)^2);
b       = -2*mu*B*(A+E);
c       = -mu^2*B^2;
fprintf('The a value is %.13f\n',a)
fprintf('The b value is %.13f\n',b)
fprintf('The c value is %.13f\n',c)

syms r
r_2vals = double(solve(r^8+a*r^6+b*r^3+c==0,r))
disp(length(r_2vals))
for i = 1:length(r_2vals)
   if real(r_2vals(i))>0
       if imag(r_2vals(i))==0
           r_2 = r_2vals(i);
       end
   end    
end
disp(r_2)
rho_2 = A+ (mu/(r_2^3))*B;
r_2_vector = R_2 +rho_2*L_2;

P = (6*(D_31*(tau_1/tau_3)+D_21*(tau/tau_3))*r_2^3)+(mu*D_31*(tau^2-tau_1^2)*(tau_1/tau_3));
Q = (6*(D_13*(tau_3/tau_1)-D_23*(tau/tau_1))*r_2^3)+(mu*D_13*(tau^2-tau_3^2)*(tau_3/tau_1));

R_1 = (6*(r_2^3))+(mu*(tau^2-tau_3^2));
R_3 = (6*(r_2^3))+(mu*(tau^2-tau_1^2));

rho_1 = (1/D_0)*((P/R_1)-D_11)
rho_3 = (1/D_0)*((Q/R_3)-D_33)


r_1_vector = R_1 +rho_1*L_1;
r_3_vector = R_3 +rho_3*L_3;
r_1_mag = norm(r_1_vector);
r_2_mag = norm(r_2_vector);
r_3_mag = norm(r_3_vector);

N = r_1_mag*(cross(r_2_vector,r_3_vector))+r_2_mag*(cross(r_3_vector,r_1_vector))+r_3_mag*(cross(r_1_vector,r_2_vector));
magN = norm(N);
D = cross(r_1_vector,r_2_vector)+cross(r_2_vector,r_3_vector)+cross(r_3_vector,r_1_vector);
magD = norm(D);

S = r_1_vector*(r_2_mag-r_3_mag)+r_2_vector*(r_3_mag-r_1_mag)+r_3_vector*(r_1_mag-r_2_mag);
v_2_vector = sqrt(mu/(magN*magD))*((cross(D,r_2_vector))/r_2_mag+S);
v_2_mag = norm(v_2_vector);

fprintf('The P value is: %13f \n',P)
fprintf('The Q value is: %13f \n',Q)
fprintf('The N value is: %13f \n',N)
fprintf('The D value is: %13f \n',D)
fprintf('The r_1_vector value is: [%13f,%13f,%13f][km] \n',r_1_vector)
fprintf('The r_3_vector value is: [%13f,%13f,%13f][km] \n',r_3_vector)
fprintf('The v_2_vector value is: [%13f,%13f,%13f][km] \n',v_2_vector)




%Start of Problem 3

%Stores initial values in the vectors
mag_r1(1) = r_1_mag;
mag_r2(1) = r_2_mag;
mag_r3(1) = r_3_mag;

mag_v2(1) = v_2_mag;

mag_rho1(1) = rho_1;
mag_rho2(1) = rho_2;
mag_rho3(1) = rho_3;

vect_r1(1,:) = r_1_vector;
vect_r2(1,:) = r_2_vector;
vect_r3(1,:) = r_3_vector;

vect_v2(1,:) = v_2_vector

error1 = 1;
error2 = 1;
error3 = 1;

i = 1;

errormin = 10^3;

syms X C S
f_1(1) = 1-((mu/(2*r_2_mag^3))*tau_1^2);
f_3(1) = 1-(mu/(2*r_2_mag^3))*tau_3^2;
g_1(1) = tau_1 -(mu/(6*r_2_mag^3))*tau_1^3;
g_3(1) = tau_3 -(mu/(6*r_2_mag^3))*tau_3^3;

%loop to converge the values of rho to createa accurate values
while error1>errormin && error2>errormin && error3>errormin
    a_invert(i) = (2/mag_r2(i))-((mag_v2(i)^2)/mu);
    mag_v2_rad(i) = dot(vect_v2(i,:),vect_r2(i,:))/mag_r2(i);
    
    C = 1/(factorial(2))- (a_invert(i)*X^2)/factorial(4) +(a_invert(i)*X^2)^2/factorial(6)...
        -(a_invert(i)*X^2)^3/factorial(8) + (a_invert(i)*X^2)^4/factorial(10) - (a_invert(i)*X^2)^5/factorial(12)...
        +(a_invert(i)*X^2)^6/factorial(14);
    S = 1/(factorial(3))- (a_invert(i)*X^2)/factorial(5) +(a_invert(i)*X^2)^2/factorial(7)...
        -(a_invert(i)*X^2)^3/factorial(9) + (a_invert(i)*X^2)^4/factorial(11) - (a_invert(i)*X^2)^5/factorial(13)...
        +(a_invert(i)*X^2)^6/factorial(15);
    %solves for x1 and x3
    X1_val(i,:) = double(solve(sqrt(mu)*tau_1 == ((mag_r2(i)*mag_v2_rad(i))/sqrt(mu))*X^2*C+(1-a_invert(i)*mag_r2(i))*X^3*S+mag_r2(i)*X,X));
    X3_val(i,:) = double(solve(sqrt(mu)*tau_3 == ((mag_r2(i)*mag_v2_rad(i))/sqrt(mu))*X^2*C+(1-a_invert(i)*mag_r2(i))*X^3*S+mag_r2(i)*X,X));
    
    %loops to make sure the x values are real and close to 0
    k = 1;
    
    for j = 1:length(X1_val(i,:))
       if X1_val(i,j) == real(X1_val(i,j))
           X_1_val(i,k) = X1_val(i,j);
           k = k +1;
       end        
    end
    
    minabs = 1000000000000;
    
    for l = 1:length(X_1_val(i,:))
        minabs1 = abs(X_1_val(i:l));
        if minabs1>0
            if minabs1<minabs
                X_1(i) = X_1_val(i,l);
                minabs = minabs1;
            end
        end
    end
    
    k = 1;
    
    for j = 1:length(X3_val(i,:))
       if X3_val(i,j) == real(X3_val(i,j))
           X_3_val(i,k) = X3_val(i,j);
           k = k+1;
       end
    end
    
    minabs = 1000000000000;
    
    for l = 1:length(X_3_val(i,:))
        minabs1 = abs(X_3_val(i,l));
        if minabs1>0
            if minabs1<minabs
                X_3(i) = X_3_val(i,l);
            end
        end
    end
    
    %Calculates the new values for the new x1 and x2 values
    C_l(i) = 1/(factorial(2))- (a_invert(i)*X_1(i)^2)/factorial(4) +(a_invert(i)*X_1(i)^2)^2/factorial(6)
        -(a_invert(i)*X_1(i)^2)^3/factorial(8) + (a_invert(i)*X_1(i)^2)^4/factorial(10) - (a_invert(i)*X_1(i)^2)^5/factorial(12)
        +(a_invert(i)*X_1(i)^2)^6/factorial(14)
    
    S_l(i) = 1/(factorial(3))- (a_invert(i)*X_1(i)^2)/factorial(5) +(a_invert(i)*X_1(i)^2)^2/factorial(7)
        -(a_invert(i)*X_1(i)^2)^3/factorial(9) + (a_invert(i)*X_1(i)^2)^4/factorial(11) - (a_invert(i)*X_1(i)^2)^5/factorial(13)
        +(a_invert(i)*X_1(i)^2)^6/factorial(15)
    
    C_3(i) = 1/(factorial(2))- (a_invert(i)*X_3(i)^2)/factorial(4) +(a_invert(i)*X_3(i)^2)^2/factorial(6)
        -(a_invert(i)*X_3(i)^2)^3/factorial(8) + (a_invert(i)*X_3(i)^2)^4/factorial(10) - (a_invert(i)*X_3(i)^2)^5/factorial(12)
        +(a_invert(i)*X_3(i)^2)^6/factorial(14)
    
    S_3(i) = 1/(factorial(3))- (a_invert(i)*X_3(i)^2)/factorial(5) +(a_invert(i)*X_3(i)^2)^2/factorial(7)
        -(a_invert(i)*X_3(i)^2)^3/factorial(9) + (a_invert(i)*X_3(i)^2)^4/factorial(11) - (a_invert(i)*X_3(i)^2)^5/factorial(13)
        +(a_invert(i)*X_3(i)^2)^6/factorial(15)
        
    f_1(i+1) = 1-((X_1(i)^2)/mag_r2(i))*C_l(i);
    g_1(i+1) = tau_1-(1/sqrt(mu))*X_1(i)^3*S_1(i);
    f_3(i+1) = 1-((X_3(i)^2)/mag_r2(i))*C_3(i);
    g_3(i+1) = tau_3-(1/sqrt(mu))*X_3(i)^3*S_3(i);
    
    f_1_bar(i) = (f_1(i)+f_1(i+1))/2; 
    f_3_bar(i) = (f_3(i)+f_3(i+1))/2;
    g_1_bar(i) = (g_1(i)+g_1(i+1))/2;
    g_3_bar(i) = (g_3(i)+g_3(i+1))/2;
    
    c1(i) = (g_3_bar(i))/(f_1_bar(i)*g_3_bar(i)-f_3_bar(i)*g_1_bar(i));
    c3(i) = (g_1_bar(i))/(f_1_bar(i)*g_3_bar(i)-f_3_bar(i)*g_1_bar(i));
    
    mag_rho1(i+1) = (1/D_0)*(-D_11+(1/cl(i))*D_21-(c3(i)/c1(i))*D_31);
    mag_rho2(i+1) = (1/D_0)*(-c1(i)*D_12+D_22-(c3(i))*D_32);
    mag_rho3(i+1) = (1/D_0)*(-(c1(i)/c3(i))*D_13+(1/c3(i))*D_23-D_33);
    
    vect_r1(i+1,:) = R_1 +mag_rho1(i+1)*L_1;
    vect_r2(i+1,:) = R_2 +mag_rho2(i+1)*L_2;
    vect_r3(i+1,:) = R_3 +mag_rho3(i+1)*L_3;
    
    mag_r1(i+1) = nomr(vect_r1(i+1,:));
    mag_r2(i+1) = nomr(vect_r2(i+1,:));
    mag_r3(i+1) = nomr(vect_r3(i+1,:));
    
    vect_v2(i+1,:) = (f_1(i+1).*vect_r3(i+1,:)-f_3(i+1).*vect_r1(i+1,:))./(f_1(i+1)*g_3(i+1)-f_3(i+1)*g_1(i+1));
    mag_v2(i+1) = norm(vect_v2(i+1,:));
    
    error1 = abs(mag_rho1(i)-mag_rho1(i+1));
    error2 = abs(mag_rho2(i)-mag_rho2(i+1));
    error3 = abs(mag_rho3(i)-mag_rho3(i+1));
    
    i = i+1;
    disp(i)
    if i>150
        error1 = 10^-9;
        error2 = 10^-9;
        error3 = 10^-9;
    end    
end

K = [0;0;1];

Hm = cross(vect_r2(i,:),vect_v2(i,:));
magHm = norm(Hm);

n = cross(K,Hm);
magn = norm(n);

evector = (1/mu)*(((mag_v2(i)^2-(mu/mag_r2(i))).*vect_r2(i,:))-((dot(vect_r2(i,:),vect_v2(i,:))).*vect_v2(i,:)));

mage = norm(evector);

p_final = (magHm^2)/mu;
a_final = p_final/(1-mage^2);
fprintf('The a of the final orbit is: %.13f \n',a_final)
r_p = p_final/(1+mage);
r_a = p_final/(1-mage);

hk = Hm(3);
ni = n(1);
nj = n(2);
ek = evector(3);
rdv = dot(vect_r2(i,:),vect_v2(i,:));

inclin = acosd(hk/magHm);

if nj>0
   Big_Omega = acosd(ni/magn);
elseif nj<0
    Big_Omega = acosd(ni/mage)+180;
else
    Big_Omega = undefined;    
end
fprintf('The longitude of assending node = %13f \n',Big_Omega)
if ek>0
Little_Omega = acosd(dot(n,evector)/(magn*mage));
elseif ek<0
    Little_Omega = acosd(dot(n,evector)/(magn*mage))+180;
else
    Little_Omega = 0;    
end
fprintf('The argument of periapse = %13f \n',Little_Omega)

if rdv>0
    nu = acosd(dot(evector,vect_r2(i,:))/(mage*mag_r2(i)));
elseif rdv<0
    nu = acosd(dot(evector,vect_r2(i,:))/(mage*mag_r2(i)))+180;
else
    nu =0;
end
fprintf('The nu = %13f \n',nu)


disp('Initial Values')
disp(' ')
fprintf('X value is %.8f km\n',X_val)
fprintf('Z value is %.8f km\n',Z)
disp(' ')
fprintf('Unit vector of rho_1 is [%.8f, %.8f, %.8f]\n',L_1)
fprintf('Unit vector of rho_2 is [%.8f, %.8f, %.8f]\n',L_2)
fprintf('Unit vector of rho_3 is [%.8f, %.8f, %.8f]\n',L_3)
disp('')
fprintf('D_0 value is %.8f\n',D_0)
fprintf('D_11 value is %.8f km\n',D_0)
fprintf('D_21 value is %.8f km\n',D_21)
fprintf('D_31 value is %.8f km\n',D_31)
fprintf('D_12 value is %.8f km\n',D_12)
fprintf('D_22 value is %.8f km\n',D_22)
fprintf('D_32 value is %.8f km\n',D_32)
fprintf('D_13 value is %.8f km\n',D_13)
fprintf('D_23 value is %.8f km\n',D_23)
fprintf('D_33 value is %.8f km\n',D_33)
disp(' ')
fprintf('A value is %.8f km\n',A)
fprintf('B value is %.8f km\n',B)
fprintf('E value is %.8f km\n',E)
disp(' ')
fprintf('a value is %.8f km^2\n',a)
fprintf('b value is %.8f km^5\n',b)
fprintf('e value is %.8f km^8\n',c)
disp(' ')
fprintf('P value is %.8f km^4\n',P(1))
fprintf('Q value is %.8f km^4\n',Q(1))
disp(' ')
fprintf('R_1 value is %.8f km^3\n',R_1(1))
fprintf('R_3 value is %.8f km^3\n',R_3(1))


for j = 1:i-1
fprintf('Iteration %g\n',j)
disp(' ')
fprintf('X_1 value %.8 km^4\n',X1_val(j,:))
fprintf('X_3 value %.8 km^4\n',X3_val(j,:))
disp(' ')
fprintf('f1 value is %.8f\n',f_1(j))
fprintf('f3 value is %.8f\n',f_3(j))
fprintf('g1 value is %.8f\n',g_1(j))
fprintf('g3 value is %.8f\n',g_3(j))
disp(' ')
fprintf('f1_bar value is %.8f\n',f_1_bar(j))
fprintf('f3_bar value is %.8f\n',f_3_bar(j))
fprintf('g1_bar value is %.8f\n',g_1_bar(j))
fprintf('g3_bar value is %.8f\n',g_3_bar(j))    
disp(' ')
fprintf('c1 value is %.8\n',c1(j))
fprintf('c3 value is %.8\n',c3(j))
disp(' ')
fprintf('rho_1 value is %.8f km\n',mag_rho1(j))
fprintf('rho_2 value is %.8f km\n',mag_rho2(j))
fprintf('rho_3 value is %.8f km\n',mag_rho3(j))
disp(' ')
fprintf('r_1 vector is [%.8f,%.8f,%.8f] [km]\n',mag_rho1(j))
fprintf('r_2 vector is [%.8f,%.8f,%.8f] [km]\n',mag_rho1(j))
fprintf('r_3 vector is [%.8f,%.8f,%.8f] [km]\n',mag_rho1(j))
disp(' ')
fprintf('r_1 value is %.8f km\n',mag_r1(j))
fprintf('r_2 value is %.8f km\n',mag_r2(j))
fprintf('r_3 value is %.8f km\n',mag_r3(j))
disp(' ')
fprintf('v_2 vector is [%.8f,%.8f,%.8f] [km/s]\n',vect_v2(j,:))
fprintf('v_2 value is %.8f km/s\n',mag_v2(j))
disp(' ')
end









