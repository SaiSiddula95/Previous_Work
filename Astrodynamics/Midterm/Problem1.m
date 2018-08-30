%Saimanoj Siddula
%AerE 451 Astrodynamics II
%Problem 1
clear,clc
format longG
mu = 3.986e5;



%The values for Orbit 1 are given with:
r_p1 = 7000; %km
r_a1 = 10000; %km
new_b0 = 90; %degrees
new_c0 = new_b0+30; %degrees
delta_new0 = new_c0-new_b0;
%Calculations for the true values of the orbit 1
a1 = (r_p1+r_a1)/2;
e1 = 1-(r_p1/a1);
p1 = a1*(1-e1^2);
fprintf('a,e,p at orbit 1 are: %.13f,%.13f km,%.13f km',a1,e1,p1)
r_b0 = p1/(1+e1*cosd(new_b0));
r_c0 = p1/(1+e1*cosd(new_c0));
energy1 = -(mu)/(2*a1);

v_b0 = sqrt(2*((-mu/(2*a1))+(mu/r_b0))); %Speed of s/c,1 at Initial
fprintf('Speed of s/c,1 at Initial %.13f\n',v_b0)
v_c0 = sqrt(2*((-mu/(2*a1))+(mu/r_c0))); %Speed of s/c,2 at Initial
v_fr = sqrt(2*((-mu/(2*a1))+(mu/r_a1))); %Speed at Final Required.
fprintf('Speed of s/c,1 at Final %.13f\n',v_fr)

%%
%Solving for the time from C at initial to 180 degres new.
E_c0 = 2*atan(sqrt((1-e1)/(1+e1))*tan((new_c0*pi/180)/2));
E_c1 = 2*atan(sqrt((1-e1)/(1+e1))*tan(pi/2)); 
delta_t = sqrt((a1^3)/mu)*((E_c1-e1*sin(E_c1))-(E_c0-e1*sin(E_c0)));    %actual time in seconds.
%%  (>^^)> <(^^)> <(^^<)

%Using Lambe't Equation:
delta_newshort = 90;
DMshort = 1; %delta_new<pi
r1 = r_b0; %Defines starting position of manuvere
r2 = r_a1; %Defines the ending position of manuver
%A=sqrt((r1*r2))*sind(delta_new0)/(sqrt(1-cosd(delta_new0)));
DM =1;
A=DM*sqrt(r1*r2*(1+cosd(delta_newshort)));
Z = 0;
C = 1/2-Z/factorial(4)+Z^2/factorial(6)-...
    Z^3/factorial(8)+Z^4/factorial(10)-Z^5/factorial(12);
S = 1/factorial(3)-Z/factorial(5)+Z^2/factorial(7)-...
    Z^3/factorial(9)+Z^4/factorial(11)-Z^5/factorial(13);

y = r1+r2-(A*(1-Z*S)/(sqrt(C)));
x = sqrt(y/C);
t = ((x^3)*S+A*sqrt(y))/sqrt(mu);
%1352
counter = 1;
while (abs(t-delta_t)>(10^-4))
    CP = (1/factorial(4))+ ((2*Z)/factorial(6))...
          -(3*(Z^2)/factorial(8))+((4*Z^3)/factorial(10));
      
    SP = (1/factorial(5))+ ((2*Z)/factorial(7))...
          -(3*(Z^2)/factorial(9))+((4*Z^3)/factorial(11));    
      
      
    dtdz = (x^3)*(SP-(3*S*CP/(2*C)))+(A/8)*(((3*S*sqrt(y))/C)+(A/x)); %Newton's Short Method to find next iteration of Z
    
    Z= Z+(delta_t-t)/dtdz;
    %fprintf('The Z through iteration %g is %.13f\n',counter,t)
    C = 1/2-Z/factorial(4)+Z^2/factorial(6)-...
        Z^3/factorial(8)+Z^4/factorial(10)-Z^5/factorial(12);
       % fprintf('The C through iteration %g is %.13f\n',counter,C)

    S = 1/factorial(3)-Z/factorial(5)+Z^2/factorial(7)-...
        Z^3/factorial(9)+Z^4/factorial(11)-Z^5/factorial(13);
        %fprintf('The S through iteration %g is %.13f\n',counter,S)
    y = r1+r2-(A*(1-Z*S)/(sqrt(C)));
    x = sqrt(y/C);
    t = ((x^3)*S+A*sqrt(y))/sqrt(mu);
   % fprintf('The time through iteration %g is %.13f\n',counter,t)
    counter = counter+1;
end

r1 = [0,r_b0,0];
r2 = [-r_a1,0,0];
f = 1- y/norm(r1);
g = A*sqrt(y/mu);
gdot = 1- y/norm(r2);
v_1 = (r2-f*r1)/g;
disp(norm(v_1))
v_2 = (gdot*r2-r1)/g;
disp(norm(v_2))

%v_b0 Speed of s/c,1 at Initial
%v_fr Speed required to be in Orbit 1 at new= 180 degrees
%disp(v_b0)
%disp(v_fr)
v_b0_vector = (sqrt(mu/p1))*[-(sind(90)),(e1+cosd(90)),0];
v_fr_vector = (sqrt(mu/p1))*[-(sind(180)),(e1+cosd(180)),0];
delta_v1 = v_1-v_b0_vector;
delta_v2 = v_2-v_fr_vector;
delta_v1_mag = norm(delta_v1);
delta_v2_mag = norm(delta_v2);
sumdeltas = delta_v1_mag +delta_v2_mag;

%using e1,mu,p1,


%Finding the values for Orbit 2:
Hm_2 = cross(r1,v_1);
e_vec = (1/mu)*(((norm(v_1)^2)-(mu/norm(r1)))*r1-(dot(r1,v_1)*v_1));
e_mag = norm(e_vec);
p = (norm(Hm_2)^2)/mu;
Em_2 = ((norm(v_1)^2)/2)-(mu/norm(r1));
a_2 = -(mu)/(2*Em_2)

fprintf('a.\n')
fprintf('The semi-major axis of orbit 2 is: %.13f[km]\n',a_2)
fprintf('b.\n')
fprintf('The eccentricity vector e of orbit 2 is: [%.13f%.13f,%.13f],[P,Q,W]\n',e_vec(1),e_vec(2),e_vec(3))
fprintf('c.\n')
fprintf('The delta_V needed at point B is: [%.13f%.13f,%.13f][km/s]\n',delta_v1(1),delta_v1(2),delta_v1(3))
fprintf('d.\n')
fprintf('The delta_V needed at point A is: [%.13f%.13f,%.13f][km/s]\n',delta_v2(1),delta_v2(2),delta_v2(3))
fprintf('e.\n')
fprintf('The total delta_V needed to complete both maneuvers: %.13f[km/s]\n',sumdeltas)





