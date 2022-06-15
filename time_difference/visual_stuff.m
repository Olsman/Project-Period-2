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

% Produce parameters for default cortex model.
clear par;
par                                         = Cullen2018CortexAxon();

% Set temperature.
par.sim.temp = 21;
    
% change simulation duration
par.sim.dt.value = 1;
par.sim.tmax.value = 10;


% change node length
par.node.geo.length.value.ref       = 0.7735;
par.node.geo.length.value.vec       = par.node.geo.length.value.ref * ones(par.geo.nnode, 1);
par.node.seg.geo.length.value.ref   = par.node.geo.length.value.ref;
par.node.seg.geo.length.value.vec   = repmat(par.node.geo.length.value.vec / par.geo.nnodeseg, 1, par.geo.nnodeseg);
par =                                 CalculateLeakConductance(par);


% change periaxonal space width and g-ratio
par.myel.geo.gratio.value.ref       = 0.6888;
par.myel.geo.gratio.value.vec_ref   = par.myel.geo.gratio.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.myel.geo.peri.value.ref         = 8.487;
par.myel.geo.peri.value.vec         = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
par.myel.geo.period.value           = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
par                                 = CalculateNumberOfMyelinLamellae(par, 'max');
 

n1 = 1;
n2 = 51;

% Set stimulus location and time
loc = [53*(n1-1)+1 53*(n2-1)+1];

for i = 1:40:4001

tim = [1 i]; 


%Run the model
Model(par, loc, tim, fullfile(saveDirectory, 'SimulationResults.mat'));

%Visualizing results
load('SimulationResults.mat')
% figure()
% mesh(MEMBRANE_POTENTIAL);

str = ['Figure_at_t_', num2str(i), '(21C).png'];

figure(i)
imagesc(MEMBRANE_POTENTIAL)
title(i)
saveas(gcf,fullfile('C:\Users\michalis\OneDrive\Υπολογιστής\PRO4003-master\21C', str))

end

%figure(2)
%for i =  1:100:length(MEMBRANE_POTENTIAL(:,1))
%     plot(MEMBRANE_POTENTIAL(i,:))
%     ylim([-90 90])
%     drawnow
%end

