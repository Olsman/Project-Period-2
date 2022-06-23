% Directory of this file
thisDirectory                               = fileparts(mfilename('fullpath'));

saveDirectory                               = fullfile(thisDirectory, 'PRO4003SimulationResults');
if ~isdir(saveDirectory)
    mkdir(saveDirectory)
end

% Regenerate MAT files for the channels.
activeChannel   = McIntyre2002SlowK;        save('SavedParameters/ActiveChannels/McIntyre2002SlowK.mat','activeChannel');
activeChannel   = McIntyre2002FastNa;       save('SavedParameters/ActiveChannels/McIntyre2002FastNa.mat','activeChannel');
activeChannel   = McIntyre2002PersistentNa; save('SavedParameters/ActiveChannels/McIntyre2002PersistentNa.mat','activeChannel');
clear activeChannel

%%

% Produce parameters for default cortex model.
clear par;
par                                         = Cullen2018CortexAxon();

% Set temperature.
par.sim.temp = 37;
    
% change simulation duration
par.sim.dt.value = 0.1;
par.sim.tmax.value = 5;

n1 = 1;
n2 = 1;

% Set stimulus location and time
loc = [53*(n1-1)+1 53*(n2-1)+1];

for i = 1

tim = [1 i]; 

%Run the model
[MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model(par, loc, tim, fullfile(saveDirectory, 'SimulationResults.mat'));

%Visualizing results
load('SimulationResults.mat')
% figure()
% mesh(MEMBRANE_POTENTIAL);

str1 = ['Plot_at_t_', num2str(i), '(37C).png'];
str2 = ['Image_at_t_', num2str(i), '(37C).png'];

figure(1)
plot(TIME_VECTOR, MEMBRANE_POTENTIAL(:,1))
title(i)
saveas(gcf,fullfile('C:\Users\michalis\OneDrive\Υπολογιστής\PRO4003-master\DAP_37', str1))

figure(2)
imagesc(MEMBRANE_POTENTIAL)
title(i)
saveas(gcf,fullfile('C:\Users\michalis\OneDrive\Υπολογιστής\PRO4003-master\DAP_37', str2))

vel1 = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);

end


%%

% Produce parameters for default cortex model.
clear par;
par                                         = Cullen2018CortexAxon();

% Set temperature.
par.sim.temp = 37;
    
% change simulation duration
par.sim.dt.value = 0.1;
par.sim.tmax.value = 5;


% change node length
par.node.geo.length.value.ref       = 0.7735;
par.node.geo.length.value.vec       = par.node.geo.length.value.ref * ones(par.geo.nnode, 1);
par.node.seg.geo.length.value.ref   = par.node.geo.length.value.ref;
par.node.seg.geo.length.value.vec   = repmat(par.node.geo.length.value.vec / par.geo.nnodeseg, 1, par.geo.nnodeseg);
par =                                 CalculateLeakConductance(par);

n1 = 1;
n2 = 1;

% Set stimulus location and time
loc = [53*(n1-1)+1 53*(n2-1)+1];

for i = 1

tim = [1 i]; 

%Run the model
[MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model(par, loc, tim, fullfile(saveDirectory, 'SimulationResults.mat'));

%Visualizing results
load('SimulationResults.mat')
% figure()
% mesh(MEMBRANE_POTENTIAL);

str1 = ['Plot_at_t_', num2str(i), '(37C).png'];
str2 = ['Image_at_t_', num2str(i), '(37C).png'];

figure(3)
plot(TIME_VECTOR, MEMBRANE_POTENTIAL(:,1))
title(i)
saveas(gcf,fullfile('C:\Users\michalis\OneDrive\Υπολογιστής\PRO4003-master\DAP_37_N', str1))

figure(4)
imagesc(MEMBRANE_POTENTIAL)
title(i)
saveas(gcf,fullfile('C:\Users\michalis\OneDrive\Υπολογιστής\PRO4003-master\DAP_37_N', str2))

vel2 = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);

end


%%

% Produce parameters for default cortex model.
clear par;
par                                         = Cullen2018CortexAxon();

% Set temperature.
par.sim.temp = 37;
    
% change simulation duration
par.sim.dt.value = 0.1;
par.sim.tmax.value = 5;


% change periaxonal space width 
par.myel.geo.peri.value.ref         = 8.487;
par.myel.geo.peri.value.vec         = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
par.myel.geo.period.value           = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
par                                 = CalculateNumberOfMyelinLamellae(par, 'max');

n1 = 1;
n2 = 1;

% Set stimulus location and time
loc = [53*(n1-1)+1 53*(n2-1)+1]; 

for i = 1
    
tim = [1 i]; 

%Run the model
[MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model(par, loc, tim, fullfile(saveDirectory, 'SimulationResults.mat'));

%Visualizing results
load('SimulationResults.mat')
% figure()
% mesh(MEMBRANE_POTENTIAL);

str1 = ['Plot_at_t_', num2str(i), '(37C).png'];
str2 = ['Image_at_t_', num2str(i), '(37C).png'];

figure(5)
plot(TIME_VECTOR, MEMBRANE_POTENTIAL(:,1))
title(i)
saveas(gcf,fullfile('C:\Users\michalis\OneDrive\Υπολογιστής\PRO4003-master\DAP_37_P', str1))

figure(6)
imagesc(MEMBRANE_POTENTIAL)
title(i)
saveas(gcf,fullfile('C:\Users\michalis\OneDrive\Υπολογιστής\PRO4003-master\DAP_37_P', str2))

vel3 = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);

end

%%

% Produce parameters for default cortex model.
clear par;
par                                         = Cullen2018CortexAxon();

% Set temperature.
par.sim.temp = 37;
    
% change simulation duration
par.sim.dt.value = 0.1;
par.sim.tmax.value = 5;

% change node length
par.node.geo.length.value.ref       = 0.7735;
par.node.geo.length.value.vec       = par.node.geo.length.value.ref * ones(par.geo.nnode, 1);
par.node.seg.geo.length.value.ref   = par.node.geo.length.value.ref;
par.node.seg.geo.length.value.vec   = repmat(par.node.geo.length.value.vec / par.geo.nnodeseg, 1, par.geo.nnodeseg);
par =                                 CalculateLeakConductance(par);

% change periaxonal space width
par.myel.geo.peri.value.ref         = 8.487;
par.myel.geo.peri.value.vec         = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
par.myel.geo.period.value           = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
par                                 = CalculateNumberOfMyelinLamellae(par, 'max');

n1 = 1;
n2 = 1;

% Set stimulus location and time
loc = [53*(n1-1)+1 53*(n2-1)+1];

for i = 1
    
tim = [1 i]; 

%Run the model
[MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model(par, loc, tim, fullfile(saveDirectory, 'SimulationResults.mat'));

%Visualizing results
load('SimulationResults.mat')
% figure()
% mesh(MEMBRANE_POTENTIAL);

str1 = ['Plot_at_t_', num2str(i), '(37C).png'];
str2 = ['Image_at_t_', num2str(i), '(37C).png'];

figure(7)
plot(TIME_VECTOR, MEMBRANE_POTENTIAL(:,1))
title(i)
saveas(gcf,fullfile('C:\Users\michalis\OneDrive\Υπολογιστής\PRO4003-master\DAP_37_P_N', str1))

figure(8)
imagesc(MEMBRANE_POTENTIAL)
title(i)
saveas(gcf,fullfile('C:\Users\michalis\OneDrive\Υπολογιστής\PRO4003-master\DAP_37_P_N', str2))

vel4 = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);

end

%%

for i =  1:100:length(MEMBRANE_POTENTIAL(:,1))
    plot(MEMBRANE_POTENTIAL(i,:));
    ylim([-90 90]);
    xlabel('Node Number');
    xlabel('Membrane Potential(mV)');
    str3 = ['Time(in ms) = ', num2str(TIME_VECTOR(i))];
    text(38, 85, str3, 'Color','red','FontSize',14);
    drawnow;
end