%Saimanoj Siddula
%AerE 451 Astrodynamics
%Homework 8
clear,clc
format long g

%{
1. Use the circular, restricted 3-body equations of motion (in 2 dimensions, z = 0) to
find a free-return solution for a trip from the Earth to the moon and back.
Make the following assumptions:

1. The lunar S/C starts in a circular LEO with an altitude of 322 km
2. Earth radius is 6378.1 km
3. The S/C LEO parking orbit and cislunar orbit are coplanar
4. The Moon has a circular orbit (e = 0) with a constant radius of 384 400 km
5. Lunar radius is 1738 km
6. The mass-distance constant µ of Earth-moon system is µ = 0.01215
7. The lunar angular velocity is ?moon = 2.64907 · 10?6
rad/s
8. Earth gravity constant is µEarth = 3.986 · 105 km3/s2
%}

%!!!!! EARTH NOT THE CENTER OF THE COORDINATE SYSTEM!!!!! MUST TAKE INTO
%ACCOUNT!!!!!!!!!!!!!!

%GIVEN VALUES
%Circular Orbit:
Rad_E = 6378.1; %km
alt_LEO = 322; %km
Rad_lunar = 1738; %km
R_lunar = 384400;%km
mew = .01215;
omega_moon = 2.64907*10^-6; %rad/sec
mew_Earth = 3.986*10^5; %km^3/s^2


R_sc = (Rad_E+alt_LEO)
gamma = 0; %degrees
Rsc_x = R_sc*cosd(gamma);
Rsc_y = R_sc*sind(gamma);

R_cmoon = -((1-mew)/R_lunar);
R_cEarth = ((mew)/R_lunar);

R_scg = R_cEarth - Rsc_x;
V_scLEO = sqrt(3.986*10^5/R_sc)
V_sc = V_scLEO/omega_moon;



options = odeset('RelTol',1e-10,'AbsTol',[1e-10,1e-10,1e-10,1e-10]);
time = 0:.01:.5;
%Initial Values
%position in x, positioan ay, velocity in x, velocity in y
guesses = [Rsc_x,Rsc_y,0,V_scLEO]./384400;
[T,Y] = ode45(@ODESolver3bd, time,guesses,options);

figure(1)
%grid on
plot(Y(:,1),Y(:,2))


%{
Enter the translunar trajectory with a single delta-V velocity change. Let the velocity
change be scalar in nature, ie, add speed in the velocity direction only, so that the
departure flight-path angle ?0 is 0. Then, figure that the acceptable departure speeds
from LEO are in a narrow range of about 10.6 to 10.98 km/s (in the Earth inertial frame).
The departure angle ?0 might be from about 75?
to 150?
. (There is no ?1 here because
this is not patched conics.) As you change these two parameters v0 and ?0, make small
and cautious variations only. Applying large changes to either parameter could cause
you to miss the correct solution without ever realizing you passed it by. Plot your results
every time you try a new trajectory to see what is happening
%}

%{
Your goal is to:
1. Achieve a peri-selenium altitude of about 250 to 400 km. Generally, free return
orbits must be a retrograde encounter.
2. Get within about 150–200 km Earth altitude on the return trip
%}



%{
Use the 3-body equations developed in class, where there is assumed to be no z motion.
Set z equal to zero everywhere in the sˆ1 and sˆ2 equations. If you use a RK routine, you
will need to put the 2nd-order DEs into first-order form. Then your initial conditions
will be the departure angle and velocity, calculated for the x and y coordinates of the
rotating frame. (My advice is to use a professionally written RKF solver.)
Keep in mind that positions and velocities in the rotating sˆ frame are not equivalent to
those in the inertial frame. You will have to derive and perform a coordinate transformation
on both position r and velocity v to change from the center-of-mass rotating sˆ frame
to the Earth-centered inertial frame, or vice-versa, where even the origin of coordinates
has moved. Be sure this is done right before attempting to find a free-return orbit.
Perform and report at least one confidence test of your 3-body code by putting a S/C
into a circular (or near circular) LEO and letting it propagate for a number of orbits (10
to 20). It should remain in a stable (and identical) orbit the entire time. You can also try
putting a S/C into the rotating frame at the lunar distance of 1 from the Earth center,
but far from the moon, say along the sˆ2 axis. The moon is far enough away that its
gravitational influence should be nearly negligible. Or, you could try the L3 Lagrangian
point on the positive sˆ1 axis. With zero velocity (in the rotating frame) the S/C should
remain in a stable orbit—if you have done everything right and are calculating with
enough precision (use at least 10 decimal digits of precision).
Good luck, and be patient. Test your program first to get confidence in its correctness.
Then, search—slowly and carefully—for the cislunar free-return orbit.

%}


%{
Present a full analysis and all your computations:
• confidence test results
• free-return initial conditions
• periselenium radius and altitude
• periselenium speed with respect to the moon
• closest approach to Earth on return.
Also plot the free-return orbit, showing its behavior in both the rotating sˆ system and
in the Earth Inertial system. Plot close-ups of the lunar pereselenium approach and the
Earth return approach.
%}




