% FUNCTION TO REPRODUCE WAVE SIMULATIONS
% FIGURE 7 of "The effect of inhibition on the existence of traveling wave
% solutions for a neural field model of human seizure termination",
% Gonzalez-Ram?rez L.R., Kramer M.A., under review, 2018.

% INHIBITORY PARAMETERS 

% ke    =   Synaptic Inhibitory Threshold, values determined by matching conditions.  
% gie   =   Inhibitory Strength, values between 0 and 1. Assuming gie=gii.
% alphai=   Inhibitory Timescale.
%           Options:    slow-acting inhibition (0.1), 
%                       inhibition acting at the same timescale as excitation (1)
%                       fast-acting inhibition (10).

ke=0.105;    
gie=0.0;
alphai=0.1;

X= wave_model_inhibition(ke, gie, alphai); 

%% PLOT EXCITATORY WAVE SOLUTION (FIGURE 7)

T0=2000;         % Time Units
x= 0:0.1:T0/100; % Time Discretization
y= 0:400:4000;   % Space Discretization

figure
imagesc(x,y,squeeze(transpose(X.u(:,500:4500))), [ -0.01 0.2]);
grid on
xlabel('Time (ms)','Fontsize',20)
ylabel('Distance (\mu m)','Fontsize',20 )
h=colorbar;
title(h,'Activity')
set(gca,'FontSize',16)

