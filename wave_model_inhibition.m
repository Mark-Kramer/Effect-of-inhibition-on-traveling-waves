% FUNCTION TO REPRODUCE WAVE SIMULATIONS
% FIGURE 7 of "The effect of inhibition on the existence of traveling wave
% solutions for a neural field model of human seizure termination",
% Gonzalez-Ram?rez L.R., Kramer M.A., under review, 2018.

% INPUT
% INHIBITORY PARAMETERS
% ke     -   Synaptic Inhibitory Threshold, values determined by matching conditions.
% gie    -   Inhibitory Strength, values between 0 and 1. Assuming gie=gii.
% alphai -   Inhibitory Timescale. Options: slow-acting inhibition (0.1),  inhibition acting at the same timescale as
%             excitation (1) and fast-acting inhibition (10).

% OUTPUT
% X.u   -   Excitatory Activity
% X.v   -   Inhibitory Activity

% BOUNDARY CONDITIONS
% We assume no-flux type boundary conditions so the actual activity grid is
% 500:N+500 units, where we have fixed N (number of spatial units) to 4 000
% (mn)


function X = wave_model_inhibition(theta, gie, alphai)

% FIXED PARAMETERS TO REPRODUCE FIGURE R7

% EXCITATORY PARAMETERS
alpha=1;   % EXCITATORY TIMESCALE
delta= alpha/10;  % ADAPTATION TIMESCALE
beta= 2.5;        % STRENGTH OF ADAPTATION

% CONNECTIVITY : The spatial integral is discretized assuming that there is a connectivity of 500 neighboring spatial units
% to each side
sigma= 600/40;  % SPACE UNITS ARE INITIALLY SCALED TO 40 SPACE UNITS AND LATER REESCALED IN THE APPROXIMATION OF THE INTEGRAL TERM
sigmaei= sigma; % We assume that sigma_ei = sigma_ee = 600 microns
sigmaie = 100/40; % We assume that sigma_ie = sigma_ii = 100 microns
sigmaii=sigmaie;

%INHIBITORY PARAMETERS
thetai=theta;  % Assuming inhibitory and excitatory synaptic thresholds are the same

% SIGMOID PARAMETER
a= 50;  % Steepness of the sigmoid function

% TIME AND SPACE UNITS
T=2000;   % Time Units - 20 ms
N=4000;    % Space Units - 4 mm
dt = 0.01;  % Time Discretization

% MODEL SIMULATION VARIABLES
u = zeros(T,N+1000);  % Excitatory Activity
v = zeros(T,N+1000);  % Inhibitory Activity
q = zeros(1,N+1000);  % Adaptation Activity
P = zeros(T,N+1000);  % Excitatory External Input
Q = zeros(T,N+1000);  % Inhibitory External Input

% THE FOLLOWING TERMS ARE USED FOR THE APPROXIMATION USED IN RUNGE KUTTA
% OF ORDER 4
uk1 = zeros(1,N+1000);
vk1 = zeros(1,N+1000);
qk1 = zeros(1,N+1000);

uk2 = zeros(1,N+1000);
vk2 = zeros(1,N+1000);
qk2 = zeros(1,N+1000);

uk3 = zeros(1,N+1000);
vk3 = zeros(1,N+1000);
qk3 = zeros(1,N+1000);

uk4 = zeros(1,N+1000);
vk4 = zeros(1,N+1000);
qk4 = zeros(1,N+1000);

ur1 = zeros(1,N+1000);
vr1 = zeros(1,N+1000);
qr1 = zeros(1,N+1000);

ur2 = zeros(1,N+1000);
vr2 = zeros(1,N+1000);
qr2 = zeros(1,N+1000);

ur3 = zeros(1,N+1000);
vr3 = zeros(1,N+1000);
qr3 = zeros(1,N+1000);

% DISCRETIZATION OF THE INTEGRAL TERM:
M= 100; % Number of number of cells per unit of space required to discretize integral
B = zeros(1, 10*M+1);
Bei = zeros(1, 10*M+1);
Bie = zeros(1, 10*M+1);
Bii = zeros(1, 10*M+1);

% Approximate integral as a SUM  u(i)*( Int (w(x-y),space(i-1),space(i+1))))
% E.g.  input from cell 500 will be,   u(500)*(Int(w((-y/4sigma), 199.8
% units ,200.2 units ) where units are measured x10^-5 m.

for j=1:5*M
    B(j) = (1/2)*(   exp ( - ( abs(2006 - (j+1)*(4)) /(40*sigma))) - exp(- ( abs( 2006 - j*(4) )/(40*sigma) ) ));
    B(10*M+2 -j) = B(j);
end
B(5*M+1)= (1-exp(-(2)/(40*sigma)));

for j=1:5*M
    Bei(j) = (1/2)*(   exp ( - ( abs(2006 - (j+1)*(4)) /(40*sigmaei))) - exp(- ( abs( 2006 - j*(4) )/(40*sigmaei) ) ));
    Bei(10*M+2 -j) = Bei(j);
    
    Bie(j) = (1/2)*(   exp ( - ( abs(2006 - (j+1)*(4)) /(40*sigmaie))) - exp(- ( abs( 2006 - j*(4) )/(40*sigmaie) ) ));
    Bie(10*M+2 -j) = gie*Bie(j);
    
    Bii(j) = (1/2)*(   exp ( - ( abs(2006 - (j+1)*(4)) /(40*sigmaii))) - exp(- ( abs( 2006 - j*(4) )/(40*sigmaii) ) ));
    Bii(10*M+2 -j) = gie*Bii(j);
end

Bei(5*M+1)= (1-exp(-(2)/(40*sigmaei)));
Bie(5*M+1)= gie*(1-exp(-(2)/(40*sigmaie)));
Bii(5*M+1)= gie*(1-exp(-(2)/(40*sigmaii)));

% EXTERNAL INPUT
p=50; % Strength of Excitatory External Input
qinh=0; % Strength of Inhibitory Excternal Input
P(50:350,501:570)= p; %Initial Excitatory Input (to initialize wave propagation)
Q(1:300,501:550)= qinh; %Initial Inhibitory Input (to)

for i=1:T-1      % time
    
    %Compute k1, first model iteration
    
    for j=501:N+500       % row cell #  ( except boundary )
        
        % input from neighboorhs
        
        C=dot(B,squeeze(u(i,j-500:j+500))) ;
        Cei=dot(Bei,squeeze(u(i,j-500:j+500))) ;
        Cie=dot(Bie,squeeze(v(i,j-500:j+500))) ;
        Cii=dot(Bii,squeeze(v(i,j-500:j+500))) ;
        
        xe=  ( C - Cie + P(i,j) ) ;
        xi=  ( Cei - Cii + Q(i,j) ) ;
        
        S= 1/(1 +exp( a*(theta-xe)));
        Si= 1/(1 +exp( a*(thetai-xi)));
        
        uk1(j)= (dt)*(-(alpha)*u(i,j) + (alpha)*S -(beta)*q(j) );
        vk1(j)= (dt)*(-(alphai)*v(i,j) + (alphai)*Si  );
        
        qk1(j)=(dt)*delta*(-q(j) + u(i,j) );
        
        ur1(j) = u(i,j) +(1/2)*uk1(j);
        vr1(j) = v(i,j) +(1/2)*vk1(j);
        qr1(j) = q(j) + (1/2)*qk1(j);
    end
    
    % Boundary Conditions
    for k=1:500
        ur1(501-k)= ur1(501+k) ;
        ur1(N+500+k)= ur1(N+500-k) ;
        vr1(501-k)= vr1(501+k) ;
        vr1(N+500+k)= vr1(N+500-k) ;
        
        qr1(501-k)= qr1(501+k) ;
        qr1(N+500+k)= qr1(N+500-k) ;
    end
    
    % Compute k2, second model iteration
    
    for l=501:N+500      % row cell #  ( except boundary )
        
        
        C= dot(B,squeeze(ur1(l-5*M:l+5*M))) ;
        Cei= dot(Bei,squeeze(ur1(l-5*M:l+5*M))) ;
        Cie= dot(Bie,squeeze(vr1(l-5*M:l+5*M))) ;
        Cii= dot(Bii,squeeze(vr1(l-5*M:l+5*M))) ;
        
        xe=  ( C - Cie + P(i,j) ) ;
        xi=  ( Cei - Cii + Q(i,j) ) ;
        
        S= 1/(1 +exp( a*(theta-xe)));
        Si= 1/(1 +exp( a*(thetai-xi)));
        
        uk2(l)= (dt)*(-(alpha)*ur1(l) +(alpha)* S -(beta)*qr1(l));
        vk2(l)= (dt)*(-(alphai)*vr1(l) +(alphai)*Si);
        qk2(l)= (dt)*delta*(-qr1(l) + ur1(l) );
        
        ur2(l) = u(i,l) +(1/2)*uk2(l);
        vr2(l) = v(i,l) +(1/2)*vk2(l);
        qr2(l) = q(l) + (1/2)*qk2(l);
        
        
    end
    
    % Boundary Conditions
    for k=1:500
        ur2(501-k)= ur2(501+k) ;
        ur2(N+500+k)= ur2(N+500-k) ;
        vr2(501-k)= vr2(501+k) ;
        vr2(N+500+k)= vr2(N+500-k) ;
        
        qr2(501-k)= qr2(501+k) ;
        qr2(N+500+k)= qr2(N+500-k) ;
    end
    
    % Compute k3, third model iteration
    
    for l=501:N+500       % row cell #  ( except boundary )
        
        C=dot(B,squeeze(ur2(l-5*M:l+5*M))) ;
        Cei=dot(Bei,squeeze(ur2(l-5*M:l+5*M))) ;
        Cie=dot(Bie,squeeze(vr2(l-5*M:l+5*M))) ;
        Cii=dot(Bii,squeeze(vr2(l-5*M:l+5*M))) ;
        
        
        xe=  ( C - Cie + P(i,j) ) ;
        xi=  ( Cei - Cii + Q(i,j) ) ;
        
        S= 1/(1 +exp( a*(theta-xe)));
        Si= 1/(1 +exp( a*(thetai-xi)));
        
        uk3(l)= (dt)*(-(alpha)*ur2(l) + (alpha)*S -(beta)*qr2(l));
        vk3(l)= (dt)*(-(alphai)*vr2(l) + (alphai)*Si);
        qk3(l)= (dt)*delta*(-qr2(l) + ur2(l) );
        
        ur3(l) = u(i,l) + uk3(l);
        vr3(l) = v(i,l) + vk3(l);
        qr3(l) = q(l) + qk3(l);
        
    end
    
    % Boundary Conditions
    for k=1:500
        ur3(501-k)= ur3(501+k) ;
        ur3(N+500+k)= ur3(N+500-k) ;
        vr3(501-k)= vr3(501+k) ;
        vr3(N+500+k)= vr3(N+500-k) ;
        
        qr3(501-k)= qr3(501+k) ;
        qr3(N+500+k)= qr3(N+500-k) ;
    end
    
    % Fourth model Iteration
    for l=501:N+500       % row cell #  ( except boundary )
        
        C= dot(B,squeeze(ur3(l-5*M:l+5*M))) ;
        Cei= dot(Bei,squeeze(ur3(l-5*M:l+5*M))) ;
        Cie= dot(Bie,squeeze(vr3(l-5*M:l+5*M))) ;
        Cii= dot(Bii,squeeze(vr3(l-5*M:l+5*M))) ;
        
        xe=  ( C - Cie + P(i,j) ) ;
        xi=  ( Cei - Cii + Q(i,j) ) ;
        
        S= 1/(1 +exp( a*(theta-xe)));
        Si= 1/(1 +exp( a*(thetai-xi)));
        
        uk4(l)= (dt)*(-(alpha)*ur3(l) + (alpha)*S -(beta)*qr3(l) );
        vk4(l)= (dt)*(-(alphai)*vr3(l) + (alphai)*Si );
        qk4(l)= (dt)*delta*(-qr3(l) + ur3(l) );
        
        u(i+1,l) = u(i,l) + (1/6)*( uk1(l) + 2*uk2(l) + 2*uk3(l) + uk4(l) );
        v(i+1,l) = v(i,l) + (1/6)*( vk1(l) + 2*vk2(l) + 2*vk3(l) + vk4(l) );
        q(l) = q(l) + (1/6)*( qk1(l) + 2*qk2(l) + 2*qk3(l) + qk4(l) );
        
    end
    
    % Boundary Conditions
    for k=1:500
        u(i+1,501-k)= u(i+1,501+k) ;
        u(i+1,N+500+k)= u(i+1,N+500-k) ;
        v(i+1,501-k)= v(i+1,501+k) ;
        v(i+1,N+500+k)= v(i+1,N+500-k) ;
        q(501-k)= q(501+k) ;
        q(N+500+k)= q(N+500-k) ;
    end
    
    % PLOT WAVE PROPAGATION
    plot(u(i,:),'r');
    hold on
    plot(v(i,:),'b');
    hold off
    xlim([1 N+500])
    ylim([-1 1.5]);
    title(['t=' num2str(i*dt)])
    drawnow;
    
end
% Return excitatory and inhibitory activity
X.u=u;
X.v=v;

end

