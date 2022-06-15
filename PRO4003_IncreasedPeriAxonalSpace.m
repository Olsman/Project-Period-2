% fileparts(mfilename(
% /Users/rosanolsmanx/Documents/MATLAB/Project Period 2/PRO4003-master

% Initiate directory for saving data.
thisDirectory   = '/Users/rosanolsmanx/Documents/MATLAB/Project Period 2/PRO4003-master/';
saveDirectory   = fullfile(thisDirectory,'Cullen2018_R0_Main');
if ~isdir(saveDirectory)
    mkdir(saveDirectory)
end

% Regenerate MAT files for the channels.
activeChannel   = McIntyre2002SlowK;        save('/Users/rosanolsmanx/Documents/MATLAB/Project Period 2/PRO4003-master/SavedParameters/ActiveChannels/McIntyre2002SlowK.mat','activeChannel');
activeChannel   = McIntyre2002FastNa;       save('/Users/rosanolsmanx/Documents/MATLAB/Project Period 2/PRO4003-master/SavedParameters/ActiveChannels/McIntyre2002FastNa.mat','activeChannel');
activeChannel   = McIntyre2002PersistentNa; save('/Users/rosanolsmanx/Documents/MATLAB/Project Period 2/PRO4003-master/SavedParameters/ActiveChannels/McIntyre2002PersistentNa.mat','activeChannel');
clear activeChannel

%%
% Initiate temperature.
temp            = [21 37];

% Run the model and calculate CV.
for k = 1:2
    
    % Produce parameters for default cortex model.
    clear par;
    par = Cullen2018CortexAxon();
    
    % Set temperature.
    par.sim.temp = temp(k);
    
    % change simulation duration
    par.sim.dt.value = 2;
    
    % Run all simulations varying periaxonal space width.
    j = 1;
    for psw = [0 20]
        par.myel.geo.peri.value.ref = psw;
        par.myel.geo.peri.value.vec = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
        par.myel.geo.period.value   = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
        par                         = CalculateNumberOfMyelinLamellae(par, 'max');
        %         par                         = UpdateInternodePeriaxonalSpaceWidth(par, par.myel.geo.peri.value.ref/2, [], [1, 2, 51, 52], 'min');
        [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model(par, fullfile(saveDirectory, ['Cullen2018Cortex_psw_' num2str(psw) '_' num2str(par.sim.temp) 'C.mat']));
        velocity.psw(j,1)           = psw;
        if k == 1
            velocity.psw(j,2)       = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
        else
            velocity.psw(j,3)       = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
        end
        j                           = j+1;
    end
    refresh;
  
    
    % Reset model and temperature.
    clear par;
    par             = Cullen2018CortexAxon();
    par.sim.temp    = temp(k);
end 

%% load the data from the two different pariaxonal spaces
psw0 = load('/Users/rosanolsmanx/Documents/MATLAB/Project Period 2/PRO4003-master/Cullen2018_R0_Main/Cullen2018Cortex_psw_0_21C.mat');
psw20 = load('/Users/rosanolsmanx/Documents/MATLAB/Project Period 2/PRO4003-master/Cullen2018_R0_Main/Cullen2018Cortex_psw_20_21C.mat');

%% look at the video - conclusion speed is slower at higher periaxonal space

figure()
for i =  1:25:length(psw0.MEMBRANE_POTENTIAL(:,1))
    plot(psw0.MEMBRANE_POTENTIAL(i,:))
    ylim([-90 90])
    drawnow
end

figure()
for i =  1:25:length(psw20.MEMBRANE_POTENTIAL(:,1))
    plot(psw20.MEMBRANE_POTENTIAL(i,:))
    ylim([-90 90])
    drawnow
end

%% can i do the same but now with two simulations at once??
% i changed the lox in the model. now look at two

% Set stimulus location
n = 51

% Run the model and calculate CV.
for k = 1:2
    
    % Produce parameters for default cortex model.
    clear par;
    par = Cullen2018CortexAxon();
    
    % Set temperature.
    par.sim.temp = temp(k);
    
    % change simulation duration
    par.sim.dt.value = 2;
    
    % Run all simulations varying periaxonal space width.
    j = 1;
    for psw = [0 20]
        par.myel.geo.peri.value.ref = psw;
        par.myel.geo.peri.value.vec = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
        par.myel.geo.period.value   = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
        par                         = CalculateNumberOfMyelinLamellae(par, 'max');
        %         par                         = UpdateInternodePeriaxonalSpaceWidth(par, par.myel.geo.peri.value.ref/2, [], [1, 2, 51, 52], 'min');
        
        % set locations
        loc = [1 53*(n-1)+1]
        
        % model
        [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model(par, loc, fullfile(saveDirectory, ['Cullen2018Cortex_pswLOC_' num2str(psw) '_' num2str(par.sim.temp) 'C.mat']));
        velocityLOC.psw(j,1)           = psw;
        if k == 1
            velocityLOC.psw(j,2)       = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
        else
            velocityLOC.psw(j,3)       = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
        end
        j                           = j+1;
    end
    refresh;
  
    
    % Reset model and temperature.
    clear par;
    par             = Cullen2018CortexAxon();
    par.sim.temp    = temp(k);
end


%% load the data from the two different pariaxonal spaces - two location stimulation
psw0LOC = load('/Users/rosanolsmanx/Documents/MATLAB/Project Period 2/PRO4003-master/Cullen2018_R0_Main/Cullen2018Cortex_pswLOC_0_21C.mat');
psw20LOC = load('/Users/rosanolsmanx/Documents/MATLAB/Project Period 2/PRO4003-master/Cullen2018_R0_Main/Cullen2018Cortex_pswLOC_20_21C.mat');
%% make figure with two simulations

figure()
for i =  1:10:length(psw0LOC.MEMBRANE_POTENTIAL(:,1))
    plot(psw0LOC.MEMBRANE_POTENTIAL(i,:))
    ylim([-90 90])
    drawnow
end

figure()
for i =  1:10:length(psw20LOC.MEMBRANE_POTENTIAL(:,1))
    plot(psw20LOC.MEMBRANE_POTENTIAL(i,:))
    ylim([-90 90])
    drawnow
end

%% now look at changing time for the last node stimulation?
% can i do the same but now with two simulations at once??
% i changed the lox in the model. now look at two

% Set stimulus location
n = 51

% Run the model and calculate CV.
for k = 1:2
    
    % Produce parameters for default cortex model.
    clear par;
    par = Cullen2018CortexAxon();
    
    % Set temperature.
    par.sim.temp = temp(k);
    
    % change simulation duration
    par.sim.dt.value = 2;
    
    % Run all simulations varying periaxonal space width.
    j = 1;
    for psw = [0 20]
        par.myel.geo.peri.value.ref = psw;
        par.myel.geo.peri.value.vec = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
        par.myel.geo.period.value   = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
        par                         = CalculateNumberOfMyelinLamellae(par, 'max');
        %         par                         = UpdateInternodePeriaxonalSpaceWidth(par, par.myel.geo.peri.value.ref/2, [], [1, 2, 51, 52], 'min');
        
        % set locations
        loc = [1 53*(n-1)+1]
        
        % set time
        tim = [1 15]
        
        % model
        [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = ModelTime(par, loc, tim, fullfile(saveDirectory, ['Cullen2018Cortex_pswLOCtime_' num2str(psw) '_' num2str(par.sim.temp) 'C.mat']));
        velocityLOCtime.psw(j,1)           = psw;
        if k == 1
            velocityLOCtime.psw(j,2)       = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
        else
            velocityLOCtime.psw(j,3)       = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
        end
        j                           = j+1;
    end
    refresh;
  
    
    % Reset model and temperature.
    clear par;
    par             = Cullen2018CortexAxon();
    par.sim.temp    = temp(k);
end

%% load the data from the two different pariaxonal spaces - two location stimulation
psw0LOCtime = load('/Users/rosanolsmanx/Documents/MATLAB/Project Period 2/PRO4003-master/Cullen2018_R0_Main/Cullen2018Cortex_pswLOCtime_0_21C.mat');
psw20LOCtime = load('/Users/rosanolsmanx/Documents/MATLAB/Project Period 2/PRO4003-master/Cullen2018_R0_Main/Cullen2018Cortex_pswLOCtime_20_21C.mat');

%% make figure with two simulations

figure()
for i =  1:10:length(psw0LOCtime.MEMBRANE_POTENTIAL(:,1))
    plot(psw0LOCtime.MEMBRANE_POTENTIAL(i,:))
    ylim([-90 90])
    drawnow
end

figure()
for i =  1:10:length(psw20LOCtime.MEMBRANE_POTENTIAL(:,1))
    plot(psw20LOCtime.MEMBRANE_POTENTIAL(i,:))
    ylim([-90 90])
    drawnow
end



%% change length of the node(?)


