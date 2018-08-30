%Saimanoj Siddula
%AerE 451 Midterm
%Problem 4
clear,clc
close all
%Given Values---------------------
r0 = [1131.34,-2282.343,6672.423]; %Initial Radius in Km
v0 = [-6.643051,4.3033,2.42879]; %Initial Velocity in Km
%r0 = [10000,0,0]; %Initial Radius in Km
%v0 = [0,9.2,0]; %Initial Velocity in Km
r0_mag = norm(r0); %km
v0_mag = norm(v0); %km/s

%Defining the Gravitational Parameter of Earth
mu = 3.986012*10^5; %km^3/sec^2

%Solving for a using a modified Vis-Viva Equation.
a= mu/((2*mu/r0_mag)-v0_mag^2); %km
%---------------------------------

%Setting the Dt parameter as the time after the initial observatin
%t = 3*60*60;
t = 40*60;
dt = t;

%Starts the initial error value at 1 to insure that the first iteration
%process begins
error = 1;

%Variable i is to count the number of iterations the program preforms and
%to use a s a storgate nu,ber of arrays
i=1;

%The first X and Z values are calculated and kept in the 1 cell of their
%respective arrays.
X(i) = (sqrt(mu)*t)/abs(a);
Z(i) = (X(i)^2)/a;
%Checked
counter = 1;
%Starting the process of finding a more accurate X and Z values by checking
%whether the error is with in the tolerance range.
while (error > 10^-8)
    C(i) = 1/2-Z(i)/factorial(4)+(Z(i)^2)/factorial(6)-(Z(i)^3)/factorial(8)+(Z(i)^4)/factorial(10)-(Z(i)^5)/factorial(12);
    S(i) = 1/factorial(3)-Z(i)/factorial(5)+(Z(i)^2)/factorial(7)-(Z(i)^3)/factorial(9)+(Z(i)^4)/factorial(11)-(Z(i)^5)/factorial(13);
    F(i) = (1-(r0_mag/a))*S(i)*X(i)^3+ (dot(r0,v0)/(sqrt(mu)))*C(i)*X(i)^2 +r0_mag*X(i)-sqrt(mu)*t;    
    F_p(i) = C(i)*X(i)^2+ (dot(r0,v0)/(sqrt(mu)))*(1-S(i)*Z(i))*X(i)+r0_mag*(1-C(i)*Z(i));  
    F_pp(i) = (1-r0_mag/a)*(1-S(i)*Z(i))*X(i)+ (dot(r0,v0)/(sqrt(mu)))*(1-C(i)*Z(i));
    
    delta(i) = 2*(4*F_p(i)^2-5*F(i)*F_pp(i))^.5;
    dx(i) = 5*F(i)/(F_p(i)+sign(F_p(i))*delta(i));
    
    X(i+1) = X(i)-dx(i);
    Z(i+1) = (X(i+1)^2)/a;
    error = abs(dx(i)^2/a);
    i = i+1;
    counter = counter+1;
end

    C(i) = 1/2-Z(i)/factorial(4)+(Z(i)^2)/factorial(6)-(Z(i)^3)/factorial(8)+(Z(i)^4)/factorial(10)-(Z(i)^5)/factorial(12);    
    S(i) = 1/factorial(3)-Z(i)/factorial(5)+(Z(i)^2)/factorial(7)-(Z(i)^3)/factorial(9)+(Z(i)^4)/factorial(11)-(Z(i)^5)/factorial(13);

f = 1-((X(i)^2)/r0_mag)*C(i);
g = t-1/sqrt(mu)*(X(i)^3)*S(i);
    fprintf('The f value at iteration is: %.13f\n',f)
    fprintf('The g value at iteration is: %.13f\n',g)

r = f*r0 + g*v0;
r_mag = norm(r);

f_dot = (sqrt(mu)/(r_mag*r0_mag))*(S(i)*Z(i)-1)*X(i);
g_dot = 1-((X(i)^2)/r_mag)*C(i);
fprintf('The fdot value at iteration is: %.13f /n',f_dot)
fprintf('The gdot value at iteration is: %.13f /n',g_dot)

v = f_dot*r0 +g_dot*v0
v_mag = norm(v);
fprintf('The V Vector is [%.13f,%.13f,%.13f][km/s] \n',v)

e_v = (1/mu).*((((v_mag^2)-(mu/r_mag)).*r)-(dot(r,v).*v));
e_mag = norm(e_v);

nu = acosd((dot(e_v,r))/(e_mag*r_mag))

Hm = cross(r,v);
p = norm(Hm)^2/mu;
Em = ((norm(v)^2)/2)-(mu/norm(r))
a = -(mu)/(2*Em)
e = sqrt(1-(p/a))


%Display the solutions to the problem
fprintf('X_%g = %.13f\n',1,X(1)) 
fprintf('Z_%g = %.13f\n',1,Z(1)) 
disp(' ')
for i=1:length(C)-1
   fprintf('Iteration: %g\n',i)
   fprintf('C_%g = %.13f\n',i,C(i))
   fprintf('S_%g = %.13f\n',i,S(i))
   fprintf('F_%g = %.13f\n',i,F(i))
   fprintf('Fp_%g = %.13f\n',i,F_p(i))
   fprintf('Fpp_%g = %.13f\n',i,F_pp(i))
   fprintf('Delta_%g = %.13f\n',i,delta(i))
   fprintf('dx_%g = %.13f\n',i,dx(i))
  %fprintf('error_%g = %.13f\n',i,error(i))
   fprintf('X_%g = %.13f\n',i,X(i+1))
   fprintf('Z_%g = %.13f\n',i,Z(i+1))
   disp(' ')   
end

fprintf('f = %0.13f\n',f)
fprintf('g = %0.13f\n',g)
fprintf('R = %.8f I +%.8f J +%.8f K km\n',r)
disp(' ')
fprintf('Magnitude of R is %.8f km\n',r_mag)
disp(' ')
fprintf('f_dot = %.13f/n',f_dot)
fprintf('g_dot = %.13f/n',g_dot)
disp(' ')
fprintf('V = %.8f I +%.8f J +%.8f K km/s\n',v)
fprintf('Magnitude of V is %.8f km/s\n',v_mag)
disp(' ')
fprintf('eccentricity of the orbit is [%.8f,%.8f,%.8f]\n',e_v)
disp(' ')
fprintf('The angle of nu is %8f [deg]\n',nu)





