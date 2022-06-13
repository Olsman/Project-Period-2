% Directory of this file
thisDirectory                               = fileparts(mfilename('fullpath'));

saveDirectory                               = fullfile(thisDirectory, 'PRO4003SimulationResults');
if ~isdir(saveDirectory)
    mkdir(saveDirectory)
end

defaultNodeLength                           = 1.50;
longNodeLength                              = 3.70;
peakNodeLength                              = 1.70;
shortNodeLength                             = 0.43;
matchCVNodeLength                           = 0.635;
gAct                                        = [30, 0.05, 0.8];
paranodePeriaxonalSpace                     = 0.012255632;
naDensityReductionPC                        = 34;
inodeLengthIncreasePC                       = 74;
PnodeLength                                 = 1.9;

% Produce parameters for default cortex model.
par                                         = Carcamo2017CortexAxon();

% Set the default node length. (This is already set in the
% Carcamo2017CortexAxon.m file, but reset here for clarity).
par.node.geo.length.value.vec(:)            = defaultNodeLength;

par.sim.dt.value = 1;

% Set stimulus location
for i = 1:10:1000
    loc = [1 87*(25-1)+1];
    par.stim.dur.value = i;
    %Run the model
    Model(par, loc, fullfile(saveDirectory, 'SimulationResults.mat'));

    %Visualizing results
    load('SimulationResults.mat')
    % figure()
    % mesh(MEMBRANE_POTENTIAL);

    figure()
    imagesc(MEMBRANE_POTENTIAL)
    colorbar
    title(strcat("Stimulus Duration: ", num2str(i), " us"))
end

% figure()
% for i =  1:50:length(MEMBRANE_POTENTIAL(:,1))
%     plot(MEMBRANE_POTENTIAL(i,:))
%     ylim([-90 90])
%     drawnow
% end


