%% Script Used for Acquisition of Experimental Data Using Cardiac Probe 
%% for Plane Wave Imaging for experiment with 3D skull segment and microbubbles
% Script written by Matthieu Toulemonde
% - Performing Plane wave Imaging with AM
% - Contrast MODE in RT to have a better image quality (Compounding with PI)
% - Start to acquire acquisition when we click on the acquisition button
% - DAS GPU Beamfomring adapted in external function
% - email: mtoulemo@ic.ac.uk / mengxing.tang@ic.ac.uk

% Test for real time saving
% Method were we save frame by frame
% How to load data will have to be changed
% * To start to copy data (DMA) in the middle of it, we have to change the
% order  of the receive transmission - acqNum to be adapted
% * Generate a final function to save all information (similar to the
% current one)



%%% Will need another script to automatically converting the saved data to
%%% Matlab format free of Verasonics information. To be done after the
%%% acquisition day
% NumBuffer = 5;
    %%% Update Resource Buffer number
    %%% Update Receive for each buffer, change the framenumber
    %%% Update Process Saving
    %%% Update Process GPU Beamforming
    
%% Include path
cd('G:\Vantage-4.4.1-2104211120');

% Verasonics 002
% cd('H:\4.4.1\Vantage-4.4.1-2104211120');

% Toulemonde
% cd('D:\Function\Verasonics\Vantage-4.4.1-2104211120');
activate; 

% Serialisation and savefast function path

% addpath(genpath('D:\Function\Toulemonde\Function'));
% addpath(genpath('D:\Function\Toulemonde\2020\04 Beaforming Last Version\Ultrasound_image_reconstruction-master/'));

% addpath(genpath('F:\Software\00 Function for All\01 Save Fast'));
% addpath(genpath('F:\Software\00 Function for All\02 Cuda Beamforming\Ultrasound_image_reconstruction-3D_diverging/'));

% %%% Verasonics 002
% addpath(genpath('H:\CustomScript\00 Function for All\01 Save Fast'));
% addpath(genpath('H:\CustomScript\00 Function for All\02 Cuda Beamforming\Ultrasound_image_reconstruction-3D_diverging/'));
% 
addpath(genpath('G:\CustomScript\Toulemonde\Function'));
addpath(genpath('G:\CustomScript\00 Function for All\01 Save Fast'));
addpath(genpath('G:\CustomScript\00 Function for All\02 Cuda Beamforming\Ultrasound_image_reconstruction-3D_diverging'));
%%
clear all;

NumBuffer = 10;
%% Saving folder
savedir=['Z:\Katharine\',datestr(now,'yymmdd'),'\'];
    if ~exist(savedir,'dir')
        mkdir(savedir);
    end
name_of_file = ['GEM5scD_Settings_' datestr(now,'yymmdd') '_' datestr(now,'HHMMSS');];

    %%% UPDATE Vantage V4
    %%% Need to initialize the variables even if modified later
    savedir_THI = savedir;
    fname_THI = [];
    
%% TGC Preset paramters
Preset.TGC(1,:)=[613,660,760,850,945,1023,1023,1023];   % does affect RF data
Preset.TGC(2,:)=[613,660,760,850,945,1023,1023,1023];   % does affect RF data

%% GPU BEAMFORMING INFORMAITON NEEDED (1/2)
UserSet.RealTime.GPUGEFull = 0;                     % If we use the full probe during Beamforming or not
UserSet.RealTime.GPUFilterChoice = 1;               % RF filter during GPU beamforming 
UserSet.RealTime.GPUFilterCoeff = [1.8 3.2; ...     % THI Frequency
                                   3.2 3.8; ...       % Optimal for BMode
                                   1.70 4.20];      % Full Probe bandwidth
UserSet.RealTime.GPURealTimeDisplayCoefficient = 5; % RealTime GPU amplification during RT
UserSet.RealTime.FrametoRecon = 1;                 % it will beamform FrametoRecon frames during GPU
% It will find the frames: ceil(linspace(1,UserSet.RealTime.numAcq,UserSet.RealTime.FrametoRecon))
    % Maximum of 10 frames of 10 angles in RT to have a 60fps
UserSet.RealTime.ReconFrameProcess = 1; %1: mean / 2:std / 3:max
UserSet.RealTime.SenseCutoff=0.3;

UserSet.RealTime.SVDEn=1;                           % SVD or Not
UserSet.RealTime.SVDEn_Random=50;                   % more 
UserSet.RealTime.ReconSkipFrame = 1; % We don't use it in realtime but maybe later with Chee Hau script

%% We just have to provide an information about the TW transmission
%%% If both are positive or one positive and negative

%%% INFORMATION
    %%% Beamforming after the acquisiion
UserSet.RealTime.PostBeamforming_Mode = 1;
    % PostBeamforming_Mode - 1: No
    % PostBeamforming_Mode - 2: Yes
    
    %%% If we start the acquisition with Bmode or PI acquisition
UserSet.RealTime.Acquisition_Mode = 1; % IN THIS SCRIPT ONLY 1
    % Acquisition_Mode - 1: BMode
    % Acquisition_Mode - 2: PI

    %%% If we beamform using BMode approach or PI approach
UserSet.RealTime.Beamforming_Mode = 1; % IN THIS SCRIPT ONLY 1
    % Beamforming_Mode - 1: BMode
    % Beamforming_Mode - 2: PI
% UserSet.RealTime.Beamforming_Mode can only be 1 when tranmission is BMode
if UserSet.RealTime.Acquisition_Mode==1
    UserSet.RealTime.Beamforming_Mode = 1;
end

%% Real-time PW acquisition information
    UserSet.RealTime.FR = 100; % Frame rate of the acquisition
    UserSet.RealTime.sampleDepthMM=50;                                
    UserSet.RealTime.TOF = ...
        ceil(UserSet.RealTime.sampleDepthMM/1000*2/1540*1e6);       % Initial TOF
    UserSet.RealTime.samplePerWave=4;    
    
%         %%% For trial - If we provide the number of frames here, it is not
%     %%% re-caculated
%         UserSet.RealTime.NumFrame = 3053; % Frames
%         UserSet.RealTime.numAcq = 1; 
%         
        %%% For trial - If we provide the number of frames here, it is not
    %%% re-caculated
        UserSet.RealTime.NumFrame = 100; % Frames
        UserSet.RealTime.numAcq = 1; 
    
    UserSet.RealTime.Sim=0;                                             % 0=hardware, 1=sim, 2=data-rerun
 
%%% Apodization for TRANSMISSION
%     apodization_FullAperture = ones(1,80); 
    apodization_FullAperture = tukeywin(80+2,0.25)'; apodization_FullAperture = apodization_FullAperture(2:end-1);
%     apodization_FullAperture = tukeywin(80+2,0.50)'; apodization_FullAperture = apodization_FullAperture(2:end-1);
UserSet.RealTime.TXApodFull = apodization_FullAperture;
    
%      apodization_HalfShape = [1 0];
   apodization_HalfShape = [1 1 0 0];
        UserSet.RealTime.TXApodHalf_value = 1;  %%% Default value
            UserSet.RealTime.TXApodHalf_All = apodization_HalfShape;
        apodization_HalfAperture = repmat(apodization_HalfShape...
            ,[1 80/size(apodization_HalfShape,2)]);
        apodization_HalfAperture = apodization_HalfAperture .* apodization_FullAperture;        
    UserSet.RealTime.TXApodHalf = apodization_HalfAperture;
%     clearvars apodization_FullAperture apodization_HalfAperture;
    
%%% Apodization for RECEPTION [during Beamforming]
%     apodization_FullAperture = ones(1,80); 
    apodization_FullAperture = tukeywin(80+2,0.25)'; apodization_FullAperture = apodization_FullAperture(2:end-1);
%     apodization_FullAperture = tukeywin(128+2,0.50)'; apodization_FullAperture = apodization_FullAperture(2:end-1);
UserSet.RealTime.RXApod = apodization_FullAperture;
%     clearvars apodization_FullAperture;
   
UserSet.RealTime.TXFreq= 2;                 % RealTime Transmit Frequency(MHz)
UserSet.RealTime.TXBW=0.67;                         % RealTime Transmit BW
UserSet.RealTime.numCycle=1;                                       % Number of half cycle

UserSet.RealTime.na = 1;            % Number of angles used 
UserSet.RealTime.AngleRange=30;     % Angle range of the transmitted angle

UserSet.RealTime.TXFocusMM=0;
UserSet.RealTime.VirtualPoint = 999;                                  % Virtual point source multiplication

UserSet.RealTime.Voltage = 8.2;     % Initial voltage when the application start

 
%%% Number of aperture during acquisition - 1/2/3 - Update 13/01/2020
%%% The position of the aperture are set automatically see Trans.HVMux.Aperture
%%% (0): 1 - Automatically use the normal aperture
%%% (1): 65 - Middle aperture of the probe
%%% (2): 1/129 - First/Last aperture of the probe. Probe split in 2
%%% (3): 1/65/129 - First/Middle/Last aperture of the probe, Verasonics example script
UserSet.RealTime.ApertureTX.value = 0;
    switch UserSet.RealTime.ApertureTX.value
        case 1
            UserSet.RealTime.ApertureTX.position = 65;
        case 2
            %%% Overlap of 16 [Loose 8 element on each side]
%             UserSet.RealTime.ApertureTX.position = [9;121];
            %%% Overlap of 32 [Loose 16 elements on each side]
            UserSet.RealTime.ApertureTX.position = [16;114];
        case 3
            %%% Verasonic default three apertures position
            UserSet.RealTime.ApertureTX.position = [1;65;129];
        otherwise
            UserSet.RealTime.ApertureTX.position = 1;
            UserSet.RealTime.ApertureTX.value = 1;
    end
    
    
%%% GE Cardiac probe has side elements that we can use or not 
UserSet.RealTime.FullArray = 1;                                    % GEM5ScD to know if we use all elements for transmission

%%% Type of acquisition, SUM (1) - PI 2 pulses / AM Full-SumHalf (2) - AM 3 pulses (3)
UserSet.RealTime.TYPE = 3; % it is used to know how to save / load data

%%% Type of Transmission, BMode (1) - PI 2 pulses (2) - AM 3 pulses (3)
UserSet.RealTime.CONTRAST = 3; %it is used for the Receive / Event
  
%% Tissue Harmonic imaging
%%% We have to initialize the GPU to use the GPU Coefficient of tissue
%%% harmonic imaging / Use different memory access (sum in the memory)
%%% Compoudning angles / PI / Different frequency
UserSet.THI = UserSet.RealTime;

UserSet.THI.GPUFilterChoice = 2;                          % We want to have by default the Second Harmonic filter when we change to THI
UserSet.THI.numAcq=1;                                          % RealTime Number of Acquisition per frame
UserSet.THI.NumFrame=4;                                        % RealTime Number of frame
    if rem(UserSet.THI.NumFrame,2)
        UserSet.THI.NumFrame=UserSet.THI.NumFrame +1;
    end
UserSet.THI.GPURealTimeDisplayCoefficient = 10;

%%% PI Tissue Harmonic Real-Time
UserSet.THI.TXFreq = [2.0];                              % RealTime BMode Transmit Frequency (MHz)
UserSet.THI.TXBW = [0.67];                               % RealTIme BMode Transmit BW

UserSet.THI.FR = 150;

UserSet.THI.numCycle=3;                                       % Number of half cycle
UserSet.THI.na = 1;
UserSet.THI.AngleRange=60;                                    % Angle range of the transmitted angle

UserSet.THI.VirtualPoint = 999;                                  % Virtual point source multiplication

%%% I have to find a way to jump from one voltage to another one when we
%%% jump from AM and THI
UserSet.THI.Voltage = 25;                                    % Initial voltage when the application start

%%% Number of aperture during acquisition - 1/2/3 - Update 13/01/2020
%%% The position of the aperture are set automatically see Trans.HVMux.Aperture
%%% (1): 65 - Middle aperture of the probe
%%% (2): 1/129 - First/Last aperture of the probe. Probe split in 2
%%% (3): 1/65/129 - First/Middle/Last aperture of the probe, Verasonics example script
UserSet.THI.ApertureTX.value = 0;
    switch UserSet.THI.ApertureTX.value
        case 1
            UserSet.THI.ApertureTX.position = 65;
        case 2
            %%% Overlap of 16 [Loose 8 element on each side]
%             UserSet.THI.ApertureTX.position = [9;121];
            %%% Overlap of 32 [Loose 16 elements on each side]
            UserSet.THI.ApertureTX.position = [16;114];
        case 3
            %%% Verasonic default three apertures position
            UserSet.THI.ApertureTX.position = [1;65;129];
        otherwise
            UserSet.THI.ApertureTX.position = 1;
            UserSet.THI.ApertureTX.value = 1;
    end

%%% GE Cardiac probe has side elements that we can use or not 
UserSet.THI.FullArray = 1;                                    % GEM5ScD to know if we use all elements for transmission

%%% Type of acquisition, SUM (1) - PI 2 pulses (2) - AM 3 pulses (3)
UserSet.THI.TYPE = 1;
    %%% WE can do the SUM or not but it is better to post process later
%%% Type of Transmission, BMode (1) - PI 2 pulses (2) - AM 3 pulses (3)
UserSet.THI.CONTRAST = 2;

%% Specify trans structure array
Trans.name = 'GEM5ScD';                                          % Probe information
Trans.units = 'mm';                                             % Unit information
Trans = computeTrans(Trans);                                    % Run Verasonics function for transducer definition
Trans.maxHighVoltage = 40;                                      % set maximum high voltage limit
Trans.AngleRange = UserSet.RealTime.AngleRange;                 % Save the AngleRange of RealTime in Trans
Trans.na = UserSet.RealTime.na;                                 % Save the number of RealTime Compounding angles in Trans

%%% WE have to update the frequency as the one provided by Verasonics in
%%% ComputeTrans is not standard. Then it mess up all memory
%%% initialization
Trans.frequency = 2.841;

%% Define system parameters
Resource.Parameters.Connector = 1;                              % Connector to use
Resource.Parameters.numTransmit = 256;                          % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;                       % number of receive channels.
Resource.Parameters.speedOfSound = 1540;                        % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.speedCorrectionFactor = 1.0;                % Speed coefficient correction
Resource.Parameters.simulateMode = UserSet.RealTime.Sim;        % Select the simulation mode
Resource.Parameters.verbose = 2;                                % warning level 0:3(error,warning,status, debug)
Resource.Parameters.initializeOnly = 0;                         % set to initialise parameter w/0 running the hardware
Resource.Parameters.simulateMode = UserSet.RealTime.Sim;

Resource.Parameters.scaleToWvl = Trans.frequency/...
    (Resource.Parameters.speedOfSound/1000);                    % MM to Wavelength

%% Definition of the transmitted angles

%%% Plane Wave transmission angles definition
if (UserSet.RealTime.na > 1)
    %%% Angles
    UserSet.RealTime.Angles = ...
        linspace(-UserSet.RealTime.AngleRange/2,...
        UserSet.RealTime.AngleRange/2,UserSet.RealTime.na)*pi/180;
    %%% Degree between each angles
    UserSet.RealTime.dtheta = UserSet.RealTime.Angles(2)-...
        UserSet.RealTime.Angles(1);
    
    %%% Triangle transmission
    temp_angles = UserSet.RealTime.Angles;          
    
    %%% In the case we have an even number of angles, we have in the
    %%% "optimal" condiction for the MoCo of J. Porée.
    if rem(UserSet.RealTime.na,2) == 0
        disp('Optimal triangle transmission');
        UserSet.RealTime.Angles(1:ceil(UserSet.RealTime.na/2)) = ...
            temp_angles(1:2:end);
        UserSet.RealTime.Angles(ceil(UserSet.RealTime.na/2)+1:end) = ...
            temp_angles(end:-2:2);
    else % Non optimal triangular transmission
        warning('Non optimal triangle transmission');
        UserSet.RealTime.Angles(1:ceil(UserSet.RealTime.na/2)) = ...
            temp_angles(1:2:end);
        UserSet.RealTime.Angles(ceil(UserSet.RealTime.na/2)+1:end) = ...
            temp_angles(end-1:-2:2);
    end
    clearvars temp_angles;
else
    UserSet.RealTime.Angles = 0;
    UserSet.RealTime.dtheta = 0;
end

%%% Diverging transmission angles for THI
if (UserSet.THI.na > 1)
    %%% Angles
    UserSet.THI.Angles = ...
        linspace(-UserSet.THI.AngleRange/2,...
        UserSet.THI.AngleRange/2,UserSet.THI.na)*pi/180;
    %%% Degree between each angles
    UserSet.THI.dtheta = UserSet.THI.Angles(2)-...
        UserSet.THI.Angles(1);
    
    %%% Triangle transmission
    temp_angles = UserSet.THI.Angles;          
    
    %%% In the case we have an even number of angles, we have in the
    %%% "optimal" condiction for the MoCo of J. Porée.
    if rem(UserSet.THI.na,2) == 0
        disp('Optimal THI triangle transmission');
        UserSet.THI.Angles(1:ceil(UserSet.THI.na/2)) = ...
            temp_angles(1:2:end);
        UserSet.THI.Angles(ceil(UserSet.THI.na/2)+1:end) = ...
            temp_angles(end:-2:2);
    else % Non optimal triangular transmission
        warning('Non optimal THI triangle transmission');
        UserSet.THI.Angles(1:ceil(UserSet.THI.na/2)) = ...
            temp_angles(1:2:end);
        UserSet.THI.Angles(ceil(UserSet.THI.na/2)+1:end) = ...
            temp_angles(end-1:-2:2);
    end
    clearvars temp_angles;
else
    UserSet.THI.Angles = 0;
    UserSet.THI.dtheta = 0;
end

% return

%% Specify TPC
TPC(1).name = 'RealTime Imaging';                          % RealTime/THI/BMode TPC
TPC(1).maxHighVoltage = 30;

%% Specify P structure array (Previous SFormat structure)
%%% RealTime and THI

%%% BEFORE AND AFTER DIVERGING and THI
%%% Allow to change the position of the virtual point source
P(1).theta = -pi/4;                                         % Angular aperture of the probe
P(1).rayDelta = 2*(-P(1).theta);                            % spacing in radians(sector) or dist. between rays
P(1).aperture = (Trans.numelements/2)*Trans.spacing;        % GEM5ScD P.aperture has to be devided by two as he only have 80 central elements
P(1).radius = (P(1).aperture/2)/tan(-P(1).theta);           % dist. to virt. apex
P(1).numRays = UserSet.RealTime.na;                         % Number of transmission 
P(1).startDepth = 0;                                        % Starting depth in wavelength
    sampleDepthWL=ceil(UserSet.RealTime.sampleDepthMM...    
        *Resource.Parameters.scaleToWvl);
P(1).endDepth = sampleDepthWL;                              % Acquisition depth in wavelengths

clearvars sampleDepthWL;


%% Maximum depth and TOF
%%% Depending of the sampling rate and maximum depth because of the
%%% diverging wave, we have to calculate again the TOF. The TOF calculate
%%% is also multiplied by a small coefficient to provide a little bit more
%%% flexibility.

% We have the maxium depth without taking into account the lateral
% propagation so we calculate the length in wavelength with it
% It is used to calculate a NEW TOF and the number of pixel we will have to
% initialized related to the maximum depth

Intermediate.maxAcqLength = sqrt(P(1).aperture^2 + ...      % Maximum distance Based on Verasonics calculation
    P(1).endDepth^2 - 2*P(1).aperture*P(1).endDepth*cos(P(1).theta-pi/2)) - P(1).startDepth;
Intermediate.wlsPer128 = Resource.Parameters.numRcvChannels/...
    (UserSet.RealTime.samplePerWave*2);                          % wavelengths in 128 samples for PI
Intermediate.numRcvSamples = 2*(Intermediate.wlsPer128*ceil(Intermediate.maxAcqLength...
    /Intermediate.wlsPer128))*UserSet.RealTime.samplePerWave;
Intermediate.maxAcqLength=Intermediate.numRcvSamples/(UserSet.RealTime.samplePerWave*2);  % final acquisition depth in wavelength


%%% ISSUE WITH THE VERASONICS, THEY CALCULATE THE NUMBER OF RCVSamples by
%%% apertureshiuft+128, so 161 elements.
% return
% We calculate the NEW TOF based on the TOF depth
NEW_TOF = ceil(Intermediate.maxAcqLength/Resource.Parameters.scaleToWvl/1000*2/1540*1e6);
% return
% If the NEW TOF is bigger, we update the new depth

%%% Update 13/01/2021, Maybe hwe have to add some extra time to the TOF to
%%% be able to change the paerture, it is coming from the Verasonics email

if NEW_TOF > UserSet.RealTime.TOF
     UserSet.RealTime.sampleDepthMM = Intermediate.maxAcqLength... % New Depth realated to the new TOF
        /Resource.Parameters.scaleToWvl;
     UserSet.THI.sampleDepthMM = UserSet.RealTime.sampleDepthMM;
     
    sampleDepthWL=ceil(UserSet.RealTime.sampleDepthMM*...  % New acquisition depth in wavelength
        Resource.Parameters.scaleToWvl);
    P(1).endDepth = sampleDepthWL*0.95;                     % PW - We multiply the acquisition depth by a coefficient to reduce RT artefacts at deep position
    P(2).endDepth = sampleDepthWL*0.95;                     % FOC - We multiply the acquisition depth by a coefficient to reduce RT artefacts at deep position
        clearvars sampleDepthWL;
        
     TOF = ceil(NEW_TOF*1.15);                              % We increase the TOF in order to be sure to have a correct time
                                                           % Change the coeffiction in order to have different TOF and reduce reflexion
        disp(['The previous TOF was : ' num2str(UserSet.RealTime.TOF) ...
            ' - New TOF is : ' num2str(TOF) ' - Depth: ' num2str(UserSet.RealTime.sampleDepthMM) ' mm']);
        
    UserSet.RealTime.TOF = TOF;                            % Save TOF in RealTime
    UserSet.THI.TOF = TOF;                            % Save TOF in THI
    
    clearvars NEW_TOF TOF;
end

%%% Update 13/01/2021 - The number of aperture used in transmission is
%%% taken into account for frame rate 

    % PRF For RealTime imaging
if (1/UserSet.RealTime.FR)*1e6 >= (UserSet.RealTime.CONTRAST*UserSet.RealTime.TOF*UserSet.RealTime.na*UserSet.RealTime.ApertureTX.value) 
    PRF_RT=(1/UserSet.RealTime.FR)*1e6-(UserSet.RealTime.CONTRAST*...
        UserSet.RealTime.ApertureTX.value*(UserSet.RealTime.na))*UserSet.RealTime.TOF;
else
    PRF_RT=0;
    UserSet.RealTime.FR= 1/(UserSet.RealTime.TOF*(UserSet.RealTime.CONTRAST*...
        UserSet.RealTime.ApertureTX.value*UserSet.RealTime.na))*1e6;
    disp(['Frame rate is reduced to ' num2str(UserSet.RealTime.FR)]);
end

    % PRF For THI Imaging
if (1/UserSet.THI.FR)*1e6 > (UserSet.THI.CONTRAST*UserSet.THI.TOF*UserSet.THI.na*UserSet.THI.ApertureTX.value) 
    PRF_THI=(1/UserSet.THI.FR)*1e6-(UserSet.THI.CONTRAST*...
        UserSet.THI.ApertureTX.value*UserSet.THI.na*UserSet.THI.TOF);
else
    PRF_THI=0;
    UserSet.THI.FR= 1/(UserSet.THI.TOF*(UserSet.THI.CONTRAST*UserSet.THI.ApertureTX.value*...
        UserSet.THI.na))*1e6;
end

% return;


%% ABout the frames, we have to calculate realated to the TOF and the frame rate 
%%% and the number of seconds. the number of frames that we need
%%% If the number of frames is already prvide, we do not caculate them
%%% again
%%% The frames for THI and the Focus have already been calculated

if ~isfield(UserSet.RealTime,'NumFrame')
    UserSet.RealTime.NumFrame = round(UserSet.RealTime.FR*UserSet.RealTime.Time/UserSet.RealTime.numAcq); %8Max
    if rem(UserSet.RealTime.NumFrame,2)
        UserSet.RealTime.NumFrame = UserSet.RealTime.NumFrame+1;
    end
end

    UserSet.RealTime.Time = UserSet.RealTime.NumFrame*UserSet.RealTime.numAcq/UserSet.RealTime.FR;
%   return  
%% Check on the GPU beamforming
if UserSet.RealTime.FrametoRecon>UserSet.RealTime.numAcq
    UserSet.RealTime.FrametoRecon = UserSet.RealTime.numAcq;
    disp(['!!! You can not process more than ' num2str(UserSet.RealTime.numAcq) ' Acquisition']);
end
%% Set up PData structure in cartesian domain
%%% Basic information
PData(1).PDelta = [0.875, 0, 0.75];
PData(1).Size(1) = 10 + ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = 10 + ceil(2*(P(1).endDepth + P(1).radius)*sin(-P(1).theta)/PData(1).PDelta(1));
PData(1).Size(3) = 1;
PData(1).Origin = [-(PData(1).Size(2)/2)*PData(1).PDelta(1),0,P(1).startDepth];

%%% RealTime / BEFORE / AFTER Diverging
PData(1).Region = struct(...
            'Shape',struct('Name','SectorFT', ...
            'Position',[0,0,-P(1).radius], ...
            'z',P(1).startDepth, ...
            'r',P(1).radius+P(1).endDepth, ...
            'angle',P(1).rayDelta, ...
            'steer',0));
PData(1).Region = computeRegions(PData(1));             % Only one region where to display 


%% Scatterers for SIMULATION
% Specify Media.  Use point targets in middle of PData(1).
% Set up Media points
% - Uncomment for speckle
% Media.MP = rand(40000,4);
% Media.MP(:,2) = 0;
% Media.MP(:,4) = 0.04*Media.MP(:,4) + 0.04;  % Random amplitude 
% Media.MP(:,1) = 2*halfwidth*(Media.MP(:,1)-0.5);
% Media.MP(:,3) = P.acqDepth*Media.MP(:,3);
pt1;
Media.attenuation = -0.5;
% Media.function = 'movePoints';

%% Ressources
%%% New method compared to previous script where we have two buffers.
%%% The first one is the acquisition (RT and acquisition)
%%% The second one is about the THI using PI transmission
%%% Update 13/01/2021 - We now can used different number of aperture in
%%% reception

%%% LOOp along buffer, need to have the right name of frames
for index_buffer = 1:NumBuffer
Resource.RcvBuffer(index_buffer).datatype = 'int16';
Resource.RcvBuffer(index_buffer).rowsPerFrame = ...                % We have all numAcq RT - WE SAVE ALL
    UserSet.RealTime.TYPE*UserSet.RealTime.ApertureTX.value*...
    Intermediate.numRcvSamples*UserSet.RealTime.numAcq*...
    UserSet.RealTime.na;
Resource.RcvBuffer(index_buffer).colsPerFrame = ...                % We have the number of channel save
    Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(index_buffer).numFrames = ...                   % We have all the NumFrame of RT
    UserSet.RealTime.NumFrame/NumBuffer;     
end
% Resources for THI imaging
Resource.RcvBuffer(NumBuffer+1).datatype = 'int16';
Resource.RcvBuffer(NumBuffer+1).rowsPerFrame = ...                % We have numAcq THI
    UserSet.THI.TYPE*UserSet.THI.ApertureTX.value*...
    Intermediate.numRcvSamples*UserSet.THI.numAcq*UserSet.THI.na;
Resource.RcvBuffer(NumBuffer+1).colsPerFrame = ...                % We have the number of channel save
    Resource.Parameters.numRcvChannels; 
Resource.RcvBuffer(NumBuffer+1).numFrames = UserSet.THI.NumFrame; % We have all the NumFrame of Destruction

%%% Definition of inter buffer
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;    % one intermediate buffer defined but not used.
Resource.InterBuffer(1).rowsPerFrame = 1024;
Resource.InterBuffer(1).colsPerFrame = PData(1).Size(2);
%%% Definition of Image Buffer
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 1;
Resource.ImageBuffer(1).rowsPerFrame = PData(1).Size(1);
Resource.ImageBuffer(1).colsPerFrame = PData(1).Size(2);

%%% In the case we load the acquistion and beamform with the Veasonics
%%% system
if UserSet.RealTime.Sim<2
    Resource.ImageBuffer(1).numFrames = 1; 
else
    Resource.ImageBuffer(1).numFrames=fix(UserSet.RealTime.numAcq*UserSet.RealTime.NumFrame/UserSet.RealTime.SimDownSample); 
end

%%% Windows Display
Resource.DisplayWindow(1).Title = 'Contrast / THI';
Resource.DisplayWindow(1).pdelta = 0.35;                 % Interpolation
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [200,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';          % Type of display - Verasonics is the new one
if UserSet.RealTime.Sim<2
    Resource.DisplayWindow(1).numFrames = 100; %Cineloop frame
else
    Resource.DisplayWindow(1).numFrames=fix(UserSet.numAcq*UserSet.NumFrame/UserSet.SimDownSample); 
end
Resource.DisplayWindow(1).AxesUnits = 'mm';
% Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).Colormap = pink(256);         % Should we change it for contrast imaging ?

%% TW structure array
% PLANE WAVE TRANSMISSION (BMode and PI)
TW(1).type = 'parametric';
TW(1).Parameters = [UserSet.RealTime.TXFreq,UserSet.RealTime.TXBW,UserSet.RealTime.numCycle,-1];   % A, B, C, D

% PLANE WAVE THI TRANSMISSION
TW(2).type = 'parametric';
TW(2).Parameters = [UserSet.THI.TXFreq(1),UserSet.THI.TXBW,UserSet.THI.numCycle,-1];   % A, B, C, D
TW(3).type = 'parametric';
TW(3).Parameters = [UserSet.THI.TXFreq(1),UserSet.THI.TXBW,UserSet.THI.numCycle,1];   % A, B, C, D
% return
  
%% Specify TX structure array.
%%% We have Coded Plane wave transmission and THI transmission
TX = repmat(struct('waveform', 1, ...                       % number of a transmit waveform structure specified above
                   'Origin', [0.0,0.0,0.0], ...             % Centre of the transducer
                   'focus', 0, ...                          % Focus point plane wave focus is at zero
                   'Steer', [0.0,0.0], ...                  % Initial steering angle
                   'Apod', zeros(1,Trans.numelements), ...  % set TX.Apod for 160 elements
                   'Delay', zeros(1,Trans.numelements),...  % Initial delay
                   'peakCutOff', 0.0,...
                   'peakBLMax', 40.50), 1, UserSet.RealTime.CONTRAST*UserSet.RealTime.na*UserSet.RealTime.ApertureTX.value+...
                    UserSet.THI.CONTRAST*UserSet.THI.na*UserSet.THI.ApertureTX.value); 
                % The number of TX depends of number of angle, contrast mode, number of aperture for RealTime / THI

%%% Update 13/01/2021
%%% How will be Transmistted the pulses
%%% Order of the transmission: Coded/Aperture/Angle...
%%% [Aperture_1_Angle_1_Half_1] [Aperture_1_Angle_1_Full][Aperture_1_Angle_1_Half_2]
%%% [Aperture_2_Angle_1_Half_1] [Aperture_2_Angle_1_Full][Aperture_2_Angle_1_Half_2]
%%% [Aperture_3_Angle_1_Half_1] [Aperture_3_Angle_1_Full][Aperture_3_Angle_1_Half_2]
%%%%% Next Angle
%%% [Aperture_1_Angle_2_Half_1] [Aperture_1_Angle_2_Full][Aperture_1_Angle_2_Half_2]
%%% [Aperture_2_Angle_2_Half_1] [Aperture_2_Angle_2_Full][Aperture_2_Angle_2_Half_2]
%%% [Aperture_3_Angle_2_Half_1] [Aperture_3_Angle_2_Full][Aperture_3_Angle_2_Half_2]
%%%%% ...
%%% RealTime

%%% Choice to have negative pulse as the MI is stronger with it
for n = 1:UserSet.RealTime.na   % na transmit events
    
    %%% Update 13/01/2021 - Loop to define all Aperture (up to 3) used if
    %%% needed.
    %%% Aperture position are defined at the beginning of the script
    for index_aperture = 1:UserSet.RealTime.ApertureTX.value        
        ktx_aperture = (index_aperture-1)*UserSet.RealTime.na;
        
        %%% First Half pulse 
        TX(n+ktx_aperture).waveform = 1;     % Corresponding TW signal
%         TX(n+ktx_aperture).aperture = UserSet.RealTime.ApertureTX.position(index_aperture); % Corresponding aperture
        TX(n+ktx_aperture).Apod(1:80) = UserSet.RealTime.TXApodHalf;
        
        TX(n+ktx_aperture).focus = -P(1).radius*UserSet.RealTime.VirtualPoint;
        
        TX(n+ktx_aperture).Steer = [UserSet.RealTime.Angles(n),0.0];    
        % We have to calculate the delays without the elements on the side
        TX(n+ktx_aperture).Delay = computeTXDelays(TX(n+ktx_aperture),'TOAE');   

        % Limit where the GE probe can start to use the side elements for GEM5ScD transmission
        if UserSet.RealTime.FullArray
            TX(n+ktx_aperture).Apod(81:160) = TX(n+ktx_aperture).Apod(1:80);
            TX(n+ktx_aperture).Delay(81:160) = TX(n+ktx_aperture).Delay(1:80);
        end
        
        %%% Second Full pulse
        k1 = UserSet.RealTime.na*UserSet.RealTime.ApertureTX.value;    
        TX(n+ktx_aperture+k1).waveform = 1;    
%         TX(n+ktx_aperture+k1).aperture = UserSet.RealTime.ApertureTX.position(index_aperture); % Corresponding aperture
        TX(n+ktx_aperture+k1).Apod(1:80) = UserSet.RealTime.TXApodFull;
        
        TX(n+ktx_aperture+k1).focus = -P(1).radius*UserSet.RealTime.VirtualPoint;

        TX(n+ktx_aperture+k1).Steer = [UserSet.RealTime.Angles(n),0.0];    
        % We have to calculate the delays without the elements on the side
        TX(n+ktx_aperture+k1).Delay = computeTXDelays(TX(n+ktx_aperture+k1),'TOAE'); 
        
        % Limit where the GE probe can start to use the side elements for GEM5ScD transmission
        if UserSet.RealTime.FullArray
            TX(n+ktx_aperture+k1).Apod(81:160) = TX(n+ktx_aperture+k1).Apod(1:80);
            TX(n+ktx_aperture+k1).Delay(81:160) = TX(n+ktx_aperture+k1).Delay(1:80);
        end
        
        %%% Third pulse Second Half
        k2 = 2*UserSet.RealTime.na*UserSet.RealTime.ApertureTX.value;    
        TX(n+ktx_aperture+k2).waveform = 1;    
%         TX(n+ktx_aperture+k2).aperture = UserSet.RealTime.ApertureTX.position(index_aperture); % Corresponding aperture
        TX(n+ktx_aperture+k2).Apod(1:80) = fliplr(UserSet.RealTime.TXApodHalf);
        
        TX(n+ktx_aperture+k2).focus = -P(1).radius*UserSet.RealTime.VirtualPoint;

        TX(n+ktx_aperture+k2).Steer = [UserSet.RealTime.Angles(n),0.0];    
        % We have to calculate the delays without the elements on the side
        TX(n+ktx_aperture+k2).Delay = computeTXDelays(TX(n+ktx_aperture+k2),'TOAE'); 
        
        % Limit where the GE probe can start to use the side elements for GEM5ScD transmission
        if UserSet.RealTime.FullArray
            TX(n+ktx_aperture+k2).Apod(81:160) = TX(n+ktx_aperture+k2).Apod(1:80);
            TX(n+ktx_aperture+k2).Delay(81:160) = TX(n+ktx_aperture+k2).Delay(1:80);
        end
    end
    clearvars index_aperture;
end
clearvars k2 k1 n;

%%% THI
k = UserSet.RealTime.CONTRAST*UserSet.RealTime.na*UserSet.RealTime.ApertureTX.value;

%%% With multiplexed probe, TOAE will mess up the GPU beamforming
%%% Negative
for n = 1:UserSet.THI.na   % na transmit events
    %%% Update 13/01/2021 - Loop to define all Aperture (up to 3) used if
    %%% needed.
    %%% Aperture position are defined at the beginning of the script
    for index_aperture = 1:UserSet.THI.ApertureTX.value    
        ktx_aperture = (index_aperture-1)*UserSet.THI.na;
        
        %%% First pulse - Negative
        TX(n+ktx_aperture+k).waveform = 2;
%         TX(n+ktx_aperture+k).aperture = UserSet.THI.ApertureTX.position(index_aperture); % Corresponding aperture
        TX(n+ktx_aperture+k).Apod(1:80) = UserSet.THI.TXApodFull;
        
        TX(n+ktx_aperture+k).focus = -P(1).radius*UserSet.THI.VirtualPoint;
        
        TX(n+ktx_aperture+k).Steer = [UserSet.THI.Angles(n),0.0];
        % We have to calculate the delays without the elements on the side
        TX(n+ktx_aperture+k).Delay = computeTXDelays(TX(n+ktx_aperture+k),'TOAE');

        % Limit where the GE probe can start to use the side elements for GEM5ScD transmission
        if UserSet.RealTime.FullArray
            TX(n+ktx_aperture+k).Apod(81:160) = TX(n+ktx_aperture+k).Apod(1:80);
            TX(n+ktx_aperture+k).Delay(81:160) = TX(n+ktx_aperture+k).Delay(1:80);
        end
        
        %%% Second pulse - Positive
        k1 = UserSet.THI.na*UserSet.THI.ApertureTX.value;    
        TX(n+ktx_aperture+k+k1).waveform = 3;
%         TX(n+ktx_aperture+k+k1).aperture = UserSet.THI.ApertureTX.position(index_aperture); % Corresponding aperture
        TX(n+ktx_aperture+k+k1).Apod(1:80) = UserSet.THI.TXApodFull;
        
        TX(n+ktx_aperture+k+k1).focus = -P(1).radius*UserSet.THI.VirtualPoint;
        
        TX(n+ktx_aperture+k+k1).Steer = [UserSet.THI.Angles(n),0.0];
        % We have to calculate the delays without the elements on the side
        TX(n+ktx_aperture+k+k1).Delay = computeTXDelays(TX(n+ktx_aperture+k+k1),'TOAE');
        
        % Limit where the GE probe can start to use the side elements for GEM5ScD transmission
        if UserSet.RealTime.FullArray
            TX(n+ktx_aperture+k+k1).Apod(81:160) = TX(n+ktx_aperture+k+k1).Apod(1:80);
            TX(n+ktx_aperture+k+k1).Delay(81:160) = TX(n+ktx_aperture+k+k1).Delay(1:80);
        end
    end
    clearvars index_aperture;
end
clearvars k k1 n;

%% Specify TGC Waveform structure.
TGC = repmat(struct('CntrlPts',[],...
    'rangeMax',[],....
    'Waveform',[]),1,size(Preset.TGC,1));

% We want two TGC because we don't want to impact the TGC of BMode mode
% to THI mode
% We need to change the TGC on the GUI when we pass from Contrast and THI
for i = 1:size(Preset.TGC,1)   
    TGC(i).CntrlPts = Preset.TGC(i,:);
    TGC(i).rangeMax = P(1).endDepth;
    TGC(i).Waveform = computeTGCWaveform(TGC(i));
end
clearvars i;

%% Specify Receive structure

%%% Receive filter that we can apply on the data directly when they are
%%% receidved. 
%%% With the GPU beamforming we don't necessery need it anymore

BPF = [];

% Specify Receive structure arrays. 
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P(1).startDepth, ...
                        'endDepth', P(1).startDepth + Intermediate.wlsPer128*ceil(Intermediate.maxAcqLength/Intermediate.wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'InputFilter',BPF, ...
                        'samplesPerWave', UserSet.RealTime.samplePerWave, ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,...
                        UserSet.RealTime.CONTRAST*UserSet.RealTime.ApertureTX.value*UserSet.RealTime.na*UserSet.RealTime.numAcq*UserSet.RealTime.NumFrame+...
                        UserSet.THI.CONTRAST*UserSet.THI.ApertureTX.value*UserSet.THI.na*UserSet.THI.numAcq*UserSet.THI.NumFrame); % Real time (BEFORE) / AFTER / Focus / THI
                    
%%% Update 13/01/2021
%%% How are saved the dataset in each Frame of RcvData
%%% Order of the saved dataset:
%%% [APERTURE 1](FirstHalfReceived / FullReceived / SecondHarlfReceived) - 
%%%     [APERTURE 2] (FirstHalfReceived / FullReceived / SecondHarlfReceived) - ...
%%% Inside [FirstHalfReceived] it will be the same as a normal script. 
%%% Angles / Acquisition
%%% [Angles_Acquisition_1][Angles_Acquisition_2]...

%%% DOing this way is easier for GPU beamforming, we directly split in "X
%%% aperture" to then split the angles / 
%%% RealTime
for index_buffer = 1:NumBuffer
    k_buffer = (index_buffer-1)*UserSet.RealTime.CONTRAST*UserSet.RealTime.na*...
            UserSet.RealTime.ApertureTX.value*UserSet.RealTime.numAcq*UserSet.RealTime.NumFrame/NumBuffer;
    for i = 1:UserSet.RealTime.NumFrame/NumBuffer  % 2 acquisitions per frame
        %%% Update 13/01/2021 - Added the Aperture number we do
%         k = (i-1)*UserSet.RealTime.CONTRAST*UserSet.RealTime.na*...
%             UserSet.RealTime.ApertureTX.value*UserSet.RealTime.numAcq;

        k = k_buffer+(i-1)*UserSet.RealTime.CONTRAST*UserSet.RealTime.na*...
            UserSet.RealTime.ApertureTX.value*UserSet.RealTime.numAcq;
        
        %%% Update 13/01/2021 - Loop along the Aperture Number
         %%% First Half receive 
         %%% Initialize the received index j
         j = 1;
         %%% Loop acquisition
         for index_numAcq = 1:UserSet.RealTime.numAcq
             %%% Loop Aperture
             for index_aperture = 1:UserSet.RealTime.ApertureTX.value
                 %%% Loop angles
                 for index_na = 1:UserSet.RealTime.na
                     %%% Buffer attribution
                     Receive(k+j).bufnum = index_buffer;
                     %%% We now apply the apodization in the Beamforming part
                     Receive(k+j).Apod(1:80)= -UserSet.RealTime.RXApod;
                     Receive(k+j).Apod(81:160)= -UserSet.RealTime.RXApod;

                     %%% Change the aperture depending of the transmission
%                      Receive(k+j).aperture = UserSet.RealTime.ApertureTX.position(index_aperture); % Corresponding aperture

                     Receive(k+j).TGC = 1;

                     Receive(k+j).callMediaFunc = 1;
                     Receive(k+j).framenum = i;
                     Receive(k+j).acqNum = j;
                     Receive(k+j).mode = 0;

                     j = j+1;
                 end
             end
         end
         clearvars j;
         clearvars index_numAcq index_aperture index_na;

         %%% Update 13/01/2021 - Loop along the Aperture Number
         %%% Full receive
         k1 = UserSet.RealTime.na*UserSet.RealTime.numAcq*UserSet.RealTime.ApertureTX.value;
         %%% Initialize the received index j
         j = 1;
         %%% Loop acquisition
         for index_numAcq = 1:UserSet.RealTime.numAcq
             %%% Loop Aperture
             for index_aperture = 1:UserSet.RealTime.ApertureTX.value
                 %%% Loop angles
                 for index_na = 1:UserSet.RealTime.na
                     %%% Buffer attribution
                     Receive(k1+k+j).bufnum = index_buffer;
                     %%% We now apply the apodization in the Beamforming part
                     Receive(k1+k+j).Apod(1:80) = UserSet.RealTime.RXApod;
                     Receive(k1+k+j).Apod(81:160) = UserSet.RealTime.RXApod;

                     %%% Change the aperture depending of the transmission
    %                  Receive(k1+k+j).aperture = UserSet.RealTime.ApertureTX.position(index_aperture); % Corresponding aperture

                     Receive(k1+k+j).TGC = 1;

                     Receive(k1+k+j).framenum = i;

                     if UserSet.RealTime.TYPE == 1 % SUM IN THE MEMORY
                        Receive(k1+k+j).acqNum = j; % Same as the previous transmissions
                        Receive(k1+k+j).mode = 1;
                     else                     % Both acquisition saved
                         Receive(k1+k+j).acqNum = k1+j; % we don't accumulate in the memory, do it in the beamforming
                        Receive(k1+k+j).mode = 0;
                     end

                     j = j+1;
                 end
             end
         end
         clearvars j;
         clearvars index_numAcq index_aperture index_na;

         %%% Update 13/01/2021 - Loop along the Aperture Number
         %%% Second Half receive 
         k2 = 2*UserSet.RealTime.na*UserSet.RealTime.numAcq*UserSet.RealTime.ApertureTX.value;
         %%% Initialize the received index j
         j = 1;
         %%% Loop acquisition
         for index_numAcq = 1:UserSet.RealTime.numAcq
             %%% Loop Aperture
             for index_aperture = 1:UserSet.RealTime.ApertureTX.value
                 %%% Loop angles
                 for index_na = 1:UserSet.RealTime.na
                     %%% Buffer attribution
                     Receive(k2+k+j).bufnum = index_buffer;
                     %%% We now apply the apodization in the Beamforming part
                     Receive(k2+k+j).Apod(1:80) = -UserSet.RealTime.RXApod;
                     Receive(k2+k+j).Apod(81:160) = -UserSet.RealTime.RXApod;
                     
                     %%% Change the aperture depending of the transmission
                     %                  Receive(k2+k+j).aperture = UserSet.RealTime.ApertureTX.position(index_aperture); % Corresponding aperture

                     Receive(k2+k+j).TGC = 1;

                     Receive(k2+k+j).framenum = i;      
                     %%% Update V04
                     if (UserSet.RealTime.TYPE == 1) || (UserSet.RealTime.TYPE == 2) % SUM IN THE MEMORY
                        Receive(k2+k+j).acqNum = j; % Same as the previous transmissions
                        Receive(k2+k+j).mode = 1;
                     else                     % Both acquisition saved
                         Receive(k2+k+j).acqNum = k2+j; % we don't accumulate in the memory, do it in the beamforming
                        Receive(k2+k+j).mode = 0;
                     end

                     j = j+1;
                 end
             end
         end
         clearvars j;
         clearvars index_numAcq index_aperture index_na;

    end
    % clearvars k k1 j numacq i;
    % return
end
 %%% Update 13/01/2021 - Added the Aperture number we do
%%% THI
k2 = UserSet.RealTime.CONTRAST*UserSet.RealTime.na*...
    UserSet.RealTime.ApertureTX.value*UserSet.RealTime.numAcq*UserSet.RealTime.NumFrame;
for i = 1:UserSet.THI.NumFrame  % 2 acquisitions per frame
    %%% Update 13/01/2021 - Added the Aperture number we do
    k = (i-1)*UserSet.THI.CONTRAST*UserSet.THI.na*...
        UserSet.THI.ApertureTX.value*UserSet.THI.numAcq;
    
    %%% Update 13/01/2021 - Loop along the Aperture Number
    %%% Positive Transmission
    %%% Initialize the received index j
    j = 1;
    %%% Loop acquisition
    for index_numAcq = 1:UserSet.THI.numAcq
         %%% Loop Aperture
         for index_aperture = 1:UserSet.THI.ApertureTX.value
             %%% Loop angles
             for index_na = 1:UserSet.THI.na
                 %%% Buffer attribution
                 Receive(k2+k+j).bufnum = NumBuffer+1;
                 %%% We now apply the apodization in the Beamforming part
                 Receive(k2+k+j).Apod(1:80) = UserSet.RealTime.RXApod;
                 Receive(k2+k+j).Apod(81:160) = UserSet.RealTime.RXApod;
                 
                 %%% Change the aperture depending of the transmission
%                  Receive(k2+k+j).aperture = UserSet.THI.ApertureTX.position(index_aperture); % Corresponding aperture
                 
                 Receive(k2+k+j).TGC = 2;
                 
                 Receive(k2+k+j).callMediaFunc = 1;
                 Receive(k2+k+j).framenum = i;
                 Receive(k2+k+j).acqNum = j;
                 Receive(k2+k+j).mode = 0;
                 
                 j = j+1;
             end
         end
     end
     clearvars j;
     clearvars index_numAcq index_aperture index_na;
     
     %%% Update 13/01/2021 - Loop along the Aperture Number
    % Negative Transmission
    %%% Initialize the received index j
    j = 1;
    %%% Loop acquisition
     k1 = UserSet.THI.na*UserSet.THI.numAcq*UserSet.THI.ApertureTX.value;
    for index_numAcq = 1:UserSet.THI.numAcq
         %%% Loop Aperture
         for index_aperture = 1:UserSet.THI.ApertureTX.value
             %%% Loop angles
             for index_na = 1:UserSet.THI.na
                 %%% Buffer attribution
                 Receive(k2+k+k1+j).bufnum = NumBuffer+1;
                 %%% We now apply the apodization in the Beamforming part
                 Receive(k2+k+k1+j).Apod(1:80) = UserSet.RealTime.RXApod;
                 Receive(k2+k+k1+j).Apod(81:160) = UserSet.RealTime.RXApod;
                 %         Receive(k2+k+j).Apod= UserSet.RealTime.RXApod;
                 
                 %%% Change the aperture depending of the transmission
%                  Receive(k2+k+k1+j).aperture = UserSet.THI.ApertureTX.position(index_aperture); % Corresponding aperture
                 
                 Receive(k2+k+k1+j).TGC = 2;
        
                 Receive(k2+k+k1+j).framenum = i;
                 
                 if UserSet.THI.TYPE == 1% SUM IN THE MEMORY
                     Receive(k2+k+k1+j).acqNum = j;  % Same as the previous transmissions
                     Receive(k2+k+k1+j).mode = 1;
                 else
                     Receive(k2+k+k1+j).acqNum = k1+j; % we don't accumulate in the memory, do it in the beamforming
                     Receive(k2+k+k1+j).mode = 0;
                 end
                 
                 j = j+1;
             end
         end
     end
     clearvars j;
     clearvars index_numAcq index_aperture index_na;     
end
clearvars k2 k k1 j i;
% return;

%% Specify Recon structure arrays.
%%% THERE IS NONE OF THEM BECAUSE WE USE CHEE HAU BEAMFORMING IN RT

%% Specify Process structure array.
Preset.compressMethod = 'log'; % 'log' or  'power'   
Preset.Compress=40;
Preset.DG=15;

pers = 0;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',Preset.DG,...            % pgain is image processing gain
                         'reject',2,...      % reject level 
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod',Preset.compressMethod,...
                         'compressFactor',Preset.Compress,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
                     
                     
% EF1 is external function for RAZ focus
Process(2).classname = 'External';
Process(2).method = 'RAZParameters';
Process(2).Parameters = {'srcbuffer','none'};      

%%% Process for Real Time GPU Beamforming
length_process = length(Process);
Process_info.GPU_init = length_process+1;
Process(length_process+1).classname = 'External';
Process(length_process+1).method = 'GPU_init';
Process(length_process+1).Parameters = {'srcbuffer','none'}; % name of buffer to process.
         
length_process = length(Process);
Process_info.GPU_Beamforming = length_process+1;
Process(length_process+1).classname = 'External';
Process(length_process+1).method = 'GPU_beamforming';
Process(length_process+1).Parameters = {'srcbuffer','receive',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,...
                         'dstbuffer','image',...
                         'dstbufnum',1,...
                         'dstframenum',-2};

length_process = length(Process);
Process_info.GPU_beamforming_THI = length_process+1;
Process(length_process+1).classname = 'External';
Process(length_process+1).method = 'GPU_beamforming_THI';
Process(length_process+1).Parameters = {'srcbuffer','receive',... % name of buffer to process.
                         'srcbufnum',NumBuffer+1,...
                         'srcframenum',-1,...
                         'dstbuffer','image',...
                         'dstbufnum',1,...
                         'dstframenum',-2};
                     
length_process = length(Process);  
Process_info.ReconRF = length_process+1;
Process(length_process+1).classname = 'External';
Process(length_process+1).method = 'ReconRF';
Process(length_process+1).Parameters = {'srcbuffer','none'};

%%% The button doing the acquisition needs to automaticillay change the name and the folder of the
%%% acquisition, want we click acquisition?

length_process = length(Process);
%%% New save function 
% Save information of the acquisition using the path generated by the
% Storage_Process_Update external function
Process_info.saveRFfastData_Information = length_process+1;
Process(length_process+1).classname = 'External';
Process(length_process+1).method = 'saveRFfastData_Information'; % calls the 'RFDataBatches' function
Process(length_process+1).Parameters ={'srcbuffer','none'}; % 
   
length_process = length(Process);
% saving_process = length_process;
Process_info.saveRFfastData_Acquisition = length_process;
for index_buffer  = 1:NumBuffer
% Acquisiton Save Data Process 
Process(length_process+index_buffer).classname = 'External';
Process(length_process+index_buffer).method = 'saveRFfastData_Acquisition'; % calls the 'RFDataBatches' function
Process(length_process+index_buffer).Parameters ={'srcbuffer','receive',...
                        'srcbufnum',index_buffer,...
                        'dstbuffer','none'};
                    
end
 %% Specify SeqControl structure arrays.
% - TOF for all acquisition
SeqControl(1).command = 'timeToNextAcq';        % time between AM / Focus pulse
SeqControl(1).argument = UserSet.RealTime.TOF;  % 20 usec
% - Frame rate for RT imaging (may only use this one)
SeqControl(2).command = 'timeToNextAcq';        % time between AM frames
SeqControl(2).argument = (PRF_RT+UserSet.RealTime.TOF);
% - Frame rate of THI Imaging
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = (PRF_THI+UserSet.THI.TOF);

SeqControl(4).command = 'returnToMatlab';       % Return to Matlab

SeqControl(5).command = 'jump';                 % jump back to nTPCRealTime
SeqControl(5).argument = 1;
SeqControl(6).command = 'jump';                 % jump back to nRealTime
SeqControl(6).argument = 1;
SeqControl(7).command = 'jump';                 % jump back to nTPCTHI
SeqControl(7).argument = 1;
SeqControl(8).command = 'jump';                 % jump back to nTHI
SeqControl(8).argument = 1;
SeqControl(9).command = 'jump';                 % jump back to nAcquisition.
SeqControl(9).argument = 1;
SeqControl(10).command = 'jump';                 % jump back to nNoBeamforming or nBeamforming .
SeqControl(10).argument = 1;

SeqControl(11).command = 'setTPCProfile'; % -- Change to Profile 1 (FlashAngle(Imaging))
SeqControl(11).condition = 'immediate';
SeqControl(11).argument = 1;

% Noop are used to synchronise hardware Software or to wait the change of 
% the TPC... PAGE 115
% set a delay between imaging and destruction (changing TPC)
% (value*200nsec; max. value is 2^25 - 1. (2^25 - 1)*200ns=6.7 sec)
% value = X(ns)/200
noopTime=1e-3; %1ms (max=104.8ms) 
SeqControl(12).command = 'noop'; 
SeqControl(12).argument = noopTime/200e-9; 
SeqControl(12).condition = 'Hw&Sw';

% Syn hardware and software
% Otherwise nothing work because the software and the hardware are not
% syncrhonise
SeqControl(13).command = 'sync';
SeqControl(13).argument = 50000000;

SeqControl(14).command = 'stop'; % Stop to freeze after AFTER Acquisition

SeqControl(15).command   = 'triggerOut';
SeqControl(15).condition = 'syncNone';      % syncNone will create trigger asap
SeqControl(15).argument  = 0;               % Trigger Delay

nsc = length(SeqControl)+1; % nsc is count of SeqControl objects

%% How many frames we have to skip in order to have a 50 fps
SKIP_FRAME_RT = round(UserSet.RealTime.FR/(30*UserSet.RealTime.numAcq));

%% To allow extra time for waiting
Resource.VDAS.dmaTimeout = 10*10000; % (ms) time software sequencer will wait for 'transferToHost'
                                    % set this time long enough to permit launching the Slave sequence
                                    % and then launching the Master sequence
%%
% Specify Event structure arrays.
n = 1; % n is count of Events

% Event(n).info = 'Initialise Arduino';
% Event(n).tx = 0;         % use 1st TX structure.
% Event(n).rcv = 0;    % use 1st Rcv structure.
% Event(n).recon = 0;      % no reconstruction.
% Event(n).process = 12;    % no processing
% Event(n).seqControl = 0; % time between syn. aper. acqs.
% n = n+1;

nTPCRealTime = n;
SeqControl(5).argument = nTPCRealTime;

Event(n).info = 'Initialise GPU';
Event(n).tx = 0;         % use 1st TX structure.
Event(n).rcv = 0;    % use 1st Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = Process_info.GPU_init;    % no processing
Event(n).seqControl = 0; % time between syn. aper. acqs.
n = n+1;

Event(n).info = 'Initialize TPC RealTime + Wait';
Event(n).tx = 0;         % use 1st TX structure.
Event(n).rcv = 0;    % use 1st Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % no processing
Event(n).seqControl = [11,12,13]; % time between syn. aper. acqs.
n = n+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Real Time
%%%%%%%%%%%%%%%%%%%%
nRealTime = n;
SeqControl(6).argument = nRealTime;

%%% For the RealTime acquisition display, we don't acquire all frames, we
%%% use only NumFrameRTDisplay
for i = 1:UserSet.RealTime.NumFrame/NumBuffer  %UserSet.RealTime.NumFrameRTDisplay / UserSet.RealTime.NumFrame 
    %%% Update 13/01/2021 - Aperture number
    k = (i-1)*(UserSet.RealTime.CONTRAST*UserSet.RealTime.na)*...
        UserSet.RealTime.ApertureTX.value*UserSet.RealTime.numAcq;
    
    %%% Update 13/01/2021 - Add the loop for aperture
    %%% Loop Acquisition    
    for h=1:UserSet.RealTime.numAcq                        
        %%% Loop Aperture
%         Event(n).info = 'Trigger Out';
%         Event(n).tx = 0;        % no TX
%         Event(n).rcv = 0;       % no Rcv
%         Event(n).recon = 0;     % no Recon
%         Event(n).process = 0;
%         Event(n).seqControl = [15]; % May not oblige to wait because samll difference of voltage
%         n = n+1;
        for index_aperture = 1:UserSet.RealTime.ApertureTX.value
            l = UserSet.RealTime.ApertureTX.value*UserSet.RealTime.na*(h-1); % We save the data at the same rcv of the first half aperture
            
            krx_aperture = (index_aperture-1)*UserSet.RealTime.na;
            k1 = UserSet.RealTime.na*UserSet.RealTime.numAcq*UserSet.RealTime.ApertureTX.value;
            k2 = 2*UserSet.RealTime.na*UserSet.RealTime.numAcq*UserSet.RealTime.ApertureTX.value;
            
            ktx_aperture = (index_aperture-1)*UserSet.RealTime.na;
            ktx = (UserSet.RealTime.na)*UserSet.RealTime.ApertureTX.value;
            ktx1 = 2*(UserSet.RealTime.na)*UserSet.RealTime.ApertureTX.value;
            
            %%% Loop angles
            for j = 1:UserSet.RealTime.na

                Event(n).info = ['RT Contrast F-Half - Frame ' num2str(i) ...
                    ' - numAcq ' num2str(h) ' - aperture ' num2str(index_aperture) ...
                    ' - na ' num2str(j)];
                Event(n).tx = j+ktx_aperture;         % use 1st TX structure.
                Event(n).rcv = l+krx_aperture+k+j;      % use 1st Rcv structure.
                Event(n).recon = 0;      % no reconstruction.
                Event(n).process = 0;    % no processing
                Event(n).seqControl = 1; % time between syn. aper. acqs.
                n = n+1;

                Event(n).info = ['RT Contrast Full']; 
                Event(n).tx = j+ktx_aperture+ktx;         % use 1st TX structure.
                Event(n).rcv = l+krx_aperture+k+k1+j;      % use 1st Rcv structure.
                Event(n).recon = 0;      % no reconstruction.
                Event(n).process = 0;    % no processing
                Event(n).seqControl = 1; % time between syn. aper. acqs.
                n = n+1;

                Event(n).info = ['RT Contrast S-Half']; 
                Event(n).tx = j+ktx_aperture+ktx1;         % use 1st TX structure.
                Event(n).rcv = l+krx_aperture+k+k2+j;      % use 1st Rcv structure.
                Event(n).recon = 0;      % no reconstruction.
                Event(n).process = 0;    % no processing
                Event(n).seqControl = 1; % time between syn. aper. acqs.
                n = n+1;
            end
        end
        Event(n-1).seqControl =[2]; % Delay between frame
        
    end
    clearvars ktx ktx1 l k1 k2 h j;
    clearvars krx_aperture ktx_aperture index_aperture;
    
    % Replace last Event's seqControl value.
    Event(n-1).seqControl = [2,nsc]; % time between frames, SeqControl struct defined below.
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1; 
   

% %       if floor(i/SKIP_FRAME_RT) == i/SKIP_FRAME_RT 
          Event(n).info = 'GPU Reconstruction';
          Event(n).tx = 0;         % no transmit
          Event(n).rcv = 0;        % no rcv
          Event(n).recon = 0;      % reconstruction
          Event(n).process = Process_info.GPU_Beamforming;    % processing
          Event(n).seqControl= 0;
          
          Event(n).seqControl=[4,nsc,nsc+1];
          
          %        Event(n).seqControl = 3;
          SeqControl(nsc).command='waitForTransferComplete';
          SeqControl(nsc+1).command='markTransferProcessed';
          SeqControl(nsc).argument=nsc-1;
          SeqControl(nsc+1).argument=nsc-1;
          nsc = nsc + 2;
          
          n = n+1;
          
          Event(n).info = 'Process Image';
          Event(n).tx = 0;         % no transmit
          Event(n).rcv = 0;        % no rcv
          Event(n).recon = 0;      % reconstruction
          Event(n).process = 1;    % processing
          Event(n).seqControl = 4;
          n=n+1;
          %   else
          %       n = n+1;
% %       end
  
end

Event(n).info = 'Jump back to RealTime event';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = [6];
n = n+1;
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  THI Real Time
%%%%%%%%%%%%%%%%%%%%
nTPCTHI = n;
SeqControl(7).argument = nTPCTHI;

Event(n).info = 'Initialise THI GPU';
Event(n).tx = 0;         % use 1st TX structure.
Event(n).rcv = 0;    % use 1st Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = Process_info.GPU_init;    % no processing
Event(n).seqControl = 0; % time between syn. aper. acqs.
n = n+1;

nTHI = n;
SeqControl(8).argument = nTHI;

%%% Update 13/01/2021 - Aperture number
k2 = UserSet.RealTime.CONTRAST*UserSet.RealTime.na*UserSet.RealTime.numAcq*...
    UserSet.RealTime.NumFrame*UserSet.RealTime.ApertureTX.value;
for i = 1:UserSet.THI.NumFrame
    %%% Update 13/01/2021 - Aperture number
    k = (i-1)*UserSet.THI.CONTRAST*UserSet.THI.numAcq*UserSet.THI.na*...
        UserSet.THI.ApertureTX.value;
    
    %%% Update 13/01/2021 - Add the loop for aperture
    %%% Loop Acquisition 
    for h=1:UserSet.THI.numAcq
        %%% Loop Aperture
        for index_aperture = 1:UserSet.THI.ApertureTX.value
            l = UserSet.THI.ApertureTX.value*UserSet.THI.na*(h-1); % We save the data at the same rcv of the first half aperture

            krx_aperture = (index_aperture-1)*UserSet.THI.na;
            k1 = UserSet.THI.na*UserSet.THI.numAcq*UserSet.THI.ApertureTX.value;

            %%% Shift of the na*contrast*aperture angle from Realtime
            ktx = (UserSet.RealTime.na*UserSet.RealTime.CONTRAST*UserSet.RealTime.ApertureTX.value);

            ktx_aperture = (index_aperture-1)*UserSet.THI.na;
            ktx1 = ktx + (UserSet.THI.na)*UserSet.THI.ApertureTX.value;

            for j = 1:UserSet.THI.na
                Event(n).info = ['RT THI Positive - Frame ' num2str(i) ...
                    ' - numAcq ' num2str(h) ' - aperture ' num2str(index_aperture) ...
                    ' - na ' num2str(j)];
                Event(n).tx = j+ktx+ktx_aperture;         % use 1st TX structure.
                Event(n).rcv = k2+krx_aperture+l+k+j;      % use 1st Rcv structure.
                Event(n).recon = 0;      % no reconstruction.
                Event(n).process = 0;    % no processing
                Event(n).seqControl = 1; % time between syn. aper. acqs.
                n = n+1;

                Event(n).info = 'RT THI Negative';
                Event(n).tx = j+ktx_aperture+ktx1;         % use 1st TX structure.
                Event(n).rcv = k2+krx_aperture+l+k+k1+j;      % use 1st Rcv structure.
                Event(n).recon = 0;      % no reconstruction.
                Event(n).process = 0;    % no processing
                Event(n).seqControl = 1; % time between syn. aper. acqs.
                n = n+1;

            end
        end
        Event(n-1).seqControl =3; % Delay between frame
    end
    clearvars ktx ktx1 l k1 h j;
    clearvars krx_aperture ktx_aperture index_aperture;
    
    % Replace last Event's seqControl value.
    Event(n-1).seqControl = [3,nsc]; % time between frames, SeqControl struct defined below.
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1; 
    
   Event(n).info = 'GPU THI Reconstruction';
   Event(n).tx = 0;         % no transmit
   Event(n).rcv = 0;        % no rcv
   Event(n).recon = 0;      % reconstruction
   Event(n).process = Process_info.GPU_beamforming_THI;    % processing
   Event(n).seqControl=[4,nsc,nsc+1];
%            if floor(i/5) == i/5  
%        Event(n).seqControl = 3;
    SeqControl(nsc).command='waitForTransferComplete';
    SeqControl(nsc+1).command='markTransferProcessed';
    SeqControl(nsc).argument=nsc-1;
    SeqControl(nsc+1).argument=nsc-1;
    nsc = nsc + 2;
%            end
   n = n+1;

   Event(n).info = 'Process Image';
   Event(n).tx = 0;         % no transmit
   Event(n).rcv = 0;        % no rcv
   Event(n).recon = 0;      % reconstruction
   Event(n).process = 1;    % processing
   Event(n).seqControl = 4;
   n=n+1;

end

Event(n).info = 'Jump back to THI event';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = 8;
n = n+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Acquisition
%%%%%%%%%%%%%%%%%%%%

%%% Destruction during acquisition BUT possibility to only do that
nAcquisition=n;
SeqControl(9).argument=nAcquisition;

%%% BMode acquisition and we acquire all the frames
%%% Initialize GPU (don't know if we jump from CONTRAST or THI
Event(n).info = 'Initialise GPU';
Event(n).tx = 0;         % use 1st TX structure.
Event(n).rcv = 0;    % use 1st Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = Process_info.GPU_init;    % no processing
Event(n).seqControl = 0; % time between syn. aper. acqs.
n = n+1;

%%%  Change the TPC
Event(n).info = 'Initializa TPC RealTime + Wait';
Event(n).tx = 0;         % use 1st TX structure.
Event(n).rcv = 0;    % use 1st Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % no processing
Event(n).seqControl = [11,12,13]; % time between syn. aper. acqs.
n = n+1;

for index_buffer = 1:NumBuffer
    k_buffer = (index_buffer-1)*UserSet.RealTime.CONTRAST*UserSet.RealTime.na*...
            UserSet.RealTime.ApertureTX.value*UserSet.RealTime.numAcq*UserSet.RealTime.NumFrame/NumBuffer;
        
    for i = 1:UserSet.RealTime.NumFrame/NumBuffer
        %%% Update 13/01/2021 - Aperture number
%         k = (i-1)*(UserSet.RealTime.CONTRAST*UserSet.RealTime.na)*...
%             UserSet.RealTime.ApertureTX.value*UserSet.RealTime.numAcq;

        k = k_buffer+(i-1)*UserSet.RealTime.CONTRAST*UserSet.RealTime.na*...
            UserSet.RealTime.ApertureTX.value*UserSet.RealTime.numAcq;
        
        %%% Update 13/01/2021 - Add the loop for aperture
        %%% Loop Acquisition    
        for h=1:UserSet.RealTime.numAcq
            %%% We only update the timeStamp for the first frame of the first
            %%% acquisition. The next one will be updated by the arduino
            %%% For trial, we will do a comparison between both at the end
    %         if h ==1 && i == 1
    %             %%% Event to have the time information of the frame*numacq
    %             Event(n).info = ['Time event - Frame ' num2str(i) ...
    %                 ' - numAcq ' num2str(h)];
    %             Event(n).tx = 0;         % use 1st TX structure.
    %             Event(n).rcv = 0;      % use 1st Rcv structure.
    %             Event(n).recon = 0;      % no reconstruction.
    %             Event(n).process = 9;    % no processing
    %             Event(n).seqControl = [4]; % time between syn. aper. acqs.
    %             n = n+1;
    %         end

            %%% Loop Aperture
            for index_aperture = 1:UserSet.RealTime.ApertureTX.value
                l = UserSet.RealTime.ApertureTX.value*UserSet.RealTime.na*(h-1); % We save the data at the same rcv of the first half aperture

                krx_aperture = (index_aperture-1)*UserSet.RealTime.na;
                k1 = UserSet.RealTime.na*UserSet.RealTime.numAcq*UserSet.RealTime.ApertureTX.value;
                k2 = 2*UserSet.RealTime.na*UserSet.RealTime.numAcq*UserSet.RealTime.ApertureTX.value;

                ktx_aperture = (index_aperture-1)*UserSet.RealTime.na;
                ktx = (UserSet.RealTime.na)*UserSet.RealTime.ApertureTX.value;
                ktx1 = 2*(UserSet.RealTime.na)*UserSet.RealTime.ApertureTX.value;

                %%% Loop angles
                for j = 1:UserSet.RealTime.na

                    Event(n).info = [num2str(index_buffer) ' - Acqui Contrast F-Half - Frame ' num2str(i) ...
                        ' - numAcq ' num2str(h) ' - aperture ' num2str(index_aperture) ...
                        ' - na ' num2str(j)];
                    Event(n).tx = j+ktx_aperture;         % use 1st TX structure.
                    Event(n).rcv = l+krx_aperture+k+j;      % use 1st Rcv structure.
                    Event(n).recon = 0;      % no reconstruction.
                    Event(n).process = 0;    % no processing
                    Event(n).seqControl = [1]; % time between syn. aper. acqs.
                    if j == 1 && index_aperture == 1 % We send a Trigger for the first angle of the first aperture
                        Event(n).seqControl = [15,1]; % time between syn. aper. acqs.
                    end
                    n = n+1;

                    Event(n).info = ['Acqui Contrast Full']; 
                    Event(n).tx = j+ktx_aperture+ktx;         % use 1st TX structure.
                    Event(n).rcv = l+krx_aperture+k+k1+j;      % use 1st Rcv structure.
                    Event(n).recon = 0;      % no reconstruction.
                    Event(n).process = 0;    % no processing
                    Event(n).seqControl = 1; % time between syn. aper. acqs.
                    n = n+1;

                    Event(n).info = ['Acqui Contrast S-Half']; 
                    Event(n).tx = j+ktx_aperture+ktx1;         % use 1st TX structure.
                    Event(n).rcv = l+krx_aperture+k+k2+j;      % use 1st Rcv structure.
                    Event(n).recon = 0;      % no reconstruction.
                    Event(n).process = 0;    % no processing
                    Event(n).seqControl = 1; % time between syn. aper. acqs.
                    n = n+1;
                end
            end
            Event(n-1).seqControl =[2]; % Delay between frame

        end
        clearvars ktx ktx1 l k1 k2 h j;
        clearvars krx_aperture ktx_aperture index_aperture;

        % Replace last Event's seqControl value.
        Event(n-1).seqControl = [2,nsc]; % time between frames, SeqControl struct defined below.
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc + 1;

%           if floor(i/SKIP_FRAME_RT) == i/SKIP_FRAME_RT 
% 
%               %%% GPU Beamforming - Remove wait for transfer          
%               Event(n).info = 'GPU Reconstruction';
%               Event(n).tx = 0;         % no transmit
%               Event(n).rcv = 0;        % no rcv
%               Event(n).recon = 0;      % reconstruction
%               Event(n).process = 5+(index_buffer-1);    % processing
%               Event(n).seqControl=[4];
% 
%               Event(n).seqControl=[4,nsc,nsc+1];
%               %            if floor(i/5) == i/5
%               %        Event(n).seqControl = 3;
%               SeqControl(nsc).command='waitForTransferComplete';
%               SeqControl(nsc+1).command='markTransferProcessed';
%               SeqControl(nsc).argument=nsc-1;
%               SeqControl(nsc+1).argument=nsc-1;
%               nsc = nsc + 2;
%               %            end
%               n = n+1;
%    
% %               n = n+1;
% 
%               Event(n).info = 'Process Image';
%               Event(n).tx = 0;         % no transmit
%               Event(n).rcv = 0;        % no rcv
%               Event(n).recon = 0;      % reconstruction
%               Event(n).process = 1;    % processing
%               Event(n).seqControl = 4;
%               n=n+1;
% 
%           end


    end
clearvars k i;

 Event(n).info = ['Save Buffer - ' num2str(index_buffer)]; %calls 'SaveRFData()'
    Event(n).tx = 0; % no transmit
    Event(n).rcv = 0; % no rcv
    Event(n).recon = 0; % no reconstruction
    Event(n).process = Process_info.saveRFfastData_Acquisition+index_buffer; % save function
    Event(n).seqControl = 4; % Return to Matlab
    n=n+1;
    SeqControl(nsc).command='waitForTransferComplete';
    SeqControl(nsc+1).command='markTransferProcessed';
    SeqControl(nsc).argument=nsc-1;
    SeqControl(nsc+1).argument=nsc-1;
    nsc = nsc + 2;

end

    Event(n).info = 'Sync';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [13];
    n = n+1;

Event(n).info = 'Save Information file'; %calls 'SaveRFData()'
Event(n).tx = 0; % no transmit
Event(n).rcv = 0; % no rcv
Event(n).recon = 0; % no reconstruction
Event(n).process = Process_info.saveRFfastData_Acquisition; % save function
Event(n).seqControl = 4; % Return to Matlab
n=n+1;

    
% By default we don't have beamforming after the acquisition
Event(n).info = 'Jump Beamforming or not';
Event(n).tx = 0; % no transmit
Event(n).rcv = 0; % no rcv
Event(n).recon = 0; % no reconstruction
Event(n).process = 0; % save function
Event(n).seqControl = 10; % Return to Matlab
n=n+1;

%%% We jump the next event if we don't do the beamforming of the
%%% acquisition
nBeamforming = n;

Event(n).info = 'Acquisition Beamforming'; %calls 'SaveRFData()'
Event(n).tx = 0; % no transmit
Event(n).rcv = 0; % no rcv
Event(n).recon = 0; % no reconstruction
Event(n).process = Process_info.ReconRF; % save function
Event(n).seqControl = [4,15]; % Return to Matlab
n=n+1;

nNoBeamforming = n;
SeqControl(10).argument=nNoBeamforming;

%%% This process has an issue if we have two angles with LRI (no trial with
%%% saving data)
Event(n).info = 'Process 2 and Sync';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 2; 
Event(n).seqControl = [15,13];
n = n+1;

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = 5;
n = n+1;

RealTimeDisplay_BModeTHI_choice = 0; % Help to know if BMode/PI or THI mode (0: BMode/PI - 1: THI)

% return
%% User specified UI Control Elements
%%% UPDATE Vantage V4.4

import vsv.seq.uicontrol.*;

UI_index = 1;
% The Cell UIInterface will be used to findobj to modify the interface
% For the sliders with need the 'Tag' (ex: 'UserA1VsSlider')
% For VsButtonGroup we need the 'Label' (Corresponds to the 'Title')
% For PushButton, it doesn't work
UIInterface_index = 1;

% - Set Voltage Acquisition Imaging
UI(UI_index).Statement ='[result,hv] = setTpcProfileHighVoltage(evalin(''base'',''UserSet.RealTime.Voltage'')'',1);';
UI(UI_index+1).Statement ='hv1Sldr = findobj(''Tag'',''hv1Sldr'');';
UI(UI_index+2).Statement ='set(hv1Sldr,''Value'',hv);';
UI(UI_index+3).Statement ='hv1Value = findobj(''Tag'',''hv1Value'');';
UI(UI_index+4).Statement ='set(hv1Value,''String'',num2str(hv,''%.1f''));';
    UI_index = UI_index+5;
        
%%% UI Interface
%%% UPDATE for Vantage V4.4
% - Range Change
MinMaxVal = [64,P(1).endDepth,P(1).endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(UI_index).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
                 'Callback',@UpdateRangeSlider);
    UI_index = UI_index+1;
    UIInterface{UIInterface_index,1} = 'Tag';UIInterface{UIInterface_index,2} = [UI(UI_index-1).Control.LocationCode UI(UI_index-1).Control.Style(3:end)];
    UIInterface{UIInterface_index,3} = ['Range (',AxesUnit,')'];UIInterface_index = UIInterface_index+1;
  
% - Choose how many frames to beamform in RT
RTbeamformingframe_txt = [1, UserSet.RealTime.numAcq, UserSet.RealTime.FrametoRecon];
UI(UI_index).Control =  VsSliderControl('LocationCode','UserB5','Label','RT - Frame',...
                  'SliderMinMaxVal',RTbeamformingframe_txt,...
                  'SliderStep',[0.1,0.2],'ValueFormat','%1.0f',...
                  'Callback',@RTBeamformingFrameSlider);
    clearvars RTbeamformingframe_txt;UI_index = UI_index+1;
   UIInterface{UIInterface_index,1} = 'Tag';UIInterface{UIInterface_index,2} = [UI(UI_index-1).Control.LocationCode UI(UI_index-1).Control.Style(3:end)];
   UIInterface{UIInterface_index,3} = ['RT - Frame'];UIInterface_index = UIInterface_index+1;
   
% - Choose which kind of process we do for RT Beamforming (mean / std / max)
RTbeamformingProcess_txt = {'MEAN','STD','SVD'};
UI(UI_index).Control = VsButtonGroupControl('LocationCode','UserB4','Title','RT - Process',...
    'NumButtons',3,'Labels',RTbeamformingProcess_txt,...
    'Callback',@RTBeamformingProcessChoice);
    clearvars RTbeamformingProcess_txt;UI_index = UI_index+1;
   UIInterface{UIInterface_index,1} = 'Title';UIInterface{UIInterface_index,2} = ['RT - Process'];UIInterface_index = UIInterface_index+1;
   
%%% UI Interface related to the display part and update in RT without
%%% anything
% - Sensitivity Cutoff
UI(UI_index).Control =  VsSliderControl('LocationCode','UserC8','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,UserSet.RealTime.SenseCutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...;
                  'Callback',@SensCutoffSlider);
    UI_index = UI_index+1;
    UIInterface{UIInterface_index,1} = 'Tag';UIInterface{UIInterface_index,2} = [UI(UI_index-1).Control.LocationCode UI(UI_index-1).Control.Style(3:end)];
    UIInterface{UIInterface_index,4} = [UI(UI_index-1).Control.LocationCode 'Edit']; % When it is a slider, add this one; 
    UIInterface{UIInterface_index,3} = ['Sens. Cutoff'];UIInterface_index = UIInterface_index+1;

% - RT GPU Display Coefficient
    GPU_Coefficient = [1, 15, UserSet.RealTime.GPURealTimeDisplayCoefficient];
UI(UI_index).Control =  VsSliderControl('LocationCode','UserC7','Label','GPU - Coeff',...
                  'SliderMinMaxVal',GPU_Coefficient,...
                  'SliderStep',[0.1,0.2],'ValueFormat','%1.0f',...
                  'Callback',@GPUDisplayCoefficientSlider);
    clearvars GPU_Coefficient;UI_index = UI_index+1;
   UIInterface{UIInterface_index,1} = 'Tag';UIInterface{UIInterface_index,2} = [UI(UI_index-1).Control.LocationCode UI(UI_index-1).Control.Style(3:end)];
   UIInterface{UIInterface_index,4} = [UI(UI_index-1).Control.LocationCode 'Edit']; % When it is a slider, add this one; 
   UIInterface{UIInterface_index,3} = ['GPU - Coeff'];UIInterface_index = UIInterface_index+1;
 
% - Choose GPU Filter coeffcieient
GPUFilter_txt = cellstr([repmat('[', [size(UserSet.RealTime.GPUFilterCoeff,1) 1])...
    num2str(UserSet.RealTime.GPUFilterCoeff(:,1),'%1.2f') ...
    repmat(' ', [size(UserSet.RealTime.GPUFilterCoeff,1) 1])...
    num2str(UserSet.RealTime.GPUFilterCoeff(:,2),'%1.2f') ...
    repmat(']', [size(UserSet.RealTime.GPUFilterCoeff,1) 1])])';
UI(UI_index).Control = VsButtonGroupControl('LocationCode','UserC5','Title','GPU - Filter MHz',...
    'NumButtons',3,'Labels',GPUFilter_txt,...
    'Callback',@GPUFIlterCoefficientChoice);
    clearvars GPUFilter_txt;UI_index = UI_index+1;
   UIInterface{UIInterface_index,1} = 'Title';UIInterface{UIInterface_index,2} = ['GPU - Filter MHz'];UIInterface_index = UIInterface_index+1;

% - BMode / THI Button
UI(UI_index).Control = VsButtonControl('LocationCode','UserB2','Label','THI',...
    'Callback',@THIButton);
    UI_index = UI_index+1;
   UIInterface{UIInterface_index,1} = 'Label';UIInterface{UIInterface_index,2} = ['THI'];UIInterface_index = UIInterface_index+1;

% - Acquisition Button
UI(UI_index).Control = VsButtonControl('LocationCode','UserC2','Label','ACQUISITION',...
    'Callback',@AcquisitionButton);
    UI_index = UI_index+1;
   UIInterface{UIInterface_index,1} = 'Label';UIInterface{UIInterface_index,2} = ['ACQUISITION'];UIInterface_index = UIInterface_index+1;

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = UserSet.RealTime.numAcq/2*UserSet.RealTime.na;

%% GPU BEAMFORMING INFORMATION NEEDED (2/2)

%%% Provide which frame to beamform with which coefficient
%%% It is also updated in Realtime when we change the number of frames to
%%% BF

%%% Check if 3D beamforming or 2D by checkking the name of the probe
if strcmpi(Trans.name(1),'M')
    sizeCoeff = [1,1,1];
else
    sizeCoeff = [1,1];
end

%%% Bmode/PI Realtime
switch UserSet.RealTime.Acquisition_Mode % Which kind of acquisition
    case 1 % BMode Transmission
        %%% Only one way to beamform, so we just have to create the
        %%% a matrix coefficient of ones having the shape of 
        %%% (numSample,numChannel,numAcq) 
        UserSet.RealTime.reconFramesCoefficient = ones([sizeCoeff UserSet.RealTime.numAcq]);

        numAcq_temp = length(UserSet.RealTime.reconFramesCoefficient);
        
    case 2 % PI Transmission
        %%% Two way to beamform:
        %%%     - Bmode (have to take care of pos/neg when we sum)
        %%%     - PI (we have less acquisition)
        switch UserSet.RealTime.Beamforming_Mode
        
            case 1 % Bmode Beamforming
                %%% Create a matrix coefficient of ones having the shape of 
                %%% (numSample,numChannel,numAcq) 
                UserSet.RealTime.reconFramesCoefficient = ones([sizeCoeff UserSet.RealTime.numAcq]);
                %%% we have a PI transmission but BMode beamforming we have to
                %%% provide a coefficient of -1 to every 2 acquisition othewise we will
                %%% sum NEG and POS acquisitions
                UserSet.RealTime.reconFramesCoefficient(1:2:end) = -1;
                numAcq_temp = length(UserSet.RealTime.reconFramesCoefficient);
                
            case 2 % PI Beamforming
                %%% we have a PI transmission and PI beamforming we have to
                %%% limit the beamforming to numAcq/2
                UserSet.RealTime.reconFramesCoefficient = ones([sizeCoeff UserSet.RealTime.numAcq/2]);
                %%% If the number of acquisition to beamform is higher than we have
                %%% for PI, we reduce the number
                if UserSet.RealTime.FrametoRecon>UserSet.RealTime.numAcq/2
                    UserSet.RealTime.FrametoRecon = UserSet.RealTime.numAcq/2;
                    disp(['You are in PI Beamfomring, maximum: ' num2str(UserSet.RealTime.numAcq/2) ' frames']);
                end
                numAcq_temp = length(UserSet.RealTime.reconFramesCoefficient);
        end        
end

% We use numAcq_temp to avoid to have too many frames for PI
UserSet.RealTime.reconFrames = ceil(linspace(1,numAcq_temp,UserSet.RealTime.FrametoRecon));
%%% As we don't beamforming everything, we just need the coefficient
%%% corresponding to the corresponding acquisition
UserSet.RealTime.reconFramesCoefficient = UserSet.RealTime.reconFramesCoefficient(UserSet.RealTime.reconFrames);
clearvars numAcq_temp;

UserSet.RealTime.scanRangle = -P(1).theta;
UserSet.RealTime.TXFocusMM = -P(1).radius*(Resource.Parameters.speedOfSound/(Trans.frequency*1e3))*UserSet.RealTime.VirtualPoint;

UserSet.THI.scanRangle = -P(1).theta;
UserSet.THI.TXFocusMM = -P(1).radius*(Resource.Parameters.speedOfSound/(Trans.frequency*1e3))*UserSet.THI.VirtualPoint;

%% EXTERNAL FUNCTION
%%% UPDATE Vantage 4.4 version
EF(1).Function = vsv.seq.function.ExFunctionDef('RAZParameters', @RAZParameters);
EF(2).Function = vsv.seq.function.ExFunctionDef('saveRFfastData_Acquisition', @saveRFfastData_Acquisition);
EF(3).Function = vsv.seq.function.ExFunctionDef('GPU_init', @GPU_init);
EF(4).Function = vsv.seq.function.ExFunctionDef('GPU_beamforming', @GPU_beamforming);
EF(5).Function = vsv.seq.function.ExFunctionDef('GPU_beamforming_THI', @GPU_beamforming_THI);
EF(6).Function = vsv.seq.function.ExFunctionDef('saveRFfastData_Information', @saveRFfastData_Information);
%%
% Save all the structures to a .mat file.
save('MatFiles/GEL3_12D_HFR_PW_AM_Aperture');
filename = 'GEL3_12D_HFR_PW_AM_Aperture.mat';
VSX;
return

%% **** Callback routines to be converted by text2cell function. ****
%%% UPDATE Vantage V4.4
% **NewUI approach**

% - Range Update
function UpdateRangeSlider(~, ~, UIValue)
%- Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');

scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P(1).endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P(1).endDepth = UIValue*scaleToWvl;
        UserSet.RealTime.sampleDepthMM = UIValue;
        UserSet.THI.sampleDepthMM = UIValue;
    else
        UserSet.RealTime.sampleDepthMM = UIValue/scaleToWvl;
        UserSet.THI.sampleDepthMM = UIValue/scaleToWvl;
    end
end
assignin('base','P',P);

PData = evalin('base','PData');

PData(1).Size(1) = ceil(P(1).endDepth)/PData.PDelta(3); 
PData(1).Size(2) = ceil(P(1).aperture/PData.PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Resource.Parameters.numRcvChannels-1)/2,0,0]; % x,y,z of upper lft crnr.

PData(1).Region = computeRegions(PData);

assignin('base','PData',PData);

evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');

%%% Certainly have to adapt it if we have several TGC
RealTimeDisplay_BModeTHI_choice = evalin('base','RealTimeDisplay_BModeTHI_choice');
TGC = evalin('base','TGC');

    if RealTimeDisplay_BModeTHI_choice %1: THI
        TGC_choice = 2;
    else
        TGC_choice = 1; %0: BMode/PI
    end

for i=length(TGC)
    TGC(i).rangeMax = P(1).endDepth;
    TGC(i).Waveform = computeTGCWaveform(TGC(i));
end
assignin('base','TGC',TGC);

% evalin('base','TGC.rangeMax = P(1).endDepth;');
% evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
evalin('base',['if VDAS==1, Result = loadTgcWaveform(' num2str(TGC_choice) '); end']);

Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');

return
end 

% - RTBeamformingFrame
function RTBeamformingFrameSlider(~, ~, UIValue)
% We choose how many frames we have to beamform in RT
% It has to take into account:
% - BMode or PI Transmission
% - BMode or PI Beamforming
% - Find the corresponding frames: ceil(linsapce(1,20,UIValue))
% keyboard

RealTimeDisplay_BModeTHI_choice = evalin('base','RealTimeDisplay_BModeTHI_choice');

% Check if we hare in BMode/PI or THI mode.
if RealTimeDisplay_BModeTHI_choice
    % We can't change the transmission kind when in THI mode
    disp('!!! You can not change this value in THI mode !!!');
    return
end

UserSet = evalin('base','UserSet');
GPUParams = evalin('base','GPUParams');
Trans = evalin('base','Trans');

% To be sure to have a round number of frames
UIValue = ceil(UIValue);

% Check if we are in BMode or PI Beamforming
if isfield(UserSet,'RealTime')
    RTbeamforming = UserSet.RealTime.Beamforming_Mode;
else
    RTbeamforming = UserSet.Beamforming_Mode;
end
% If we are in PI Beamforming, the numAcq is limited to half of UserSet.RealTime.numAcq
if RTbeamforming == 2 && UIValue>GPUParams.numAcq
    UIValue = GPUParams.numAcq/2;
    disp(['You are in PI Beamfomring, maximum: ' num2str(GPUParams.numAcq) ' frames']);
end
    
if isfield(UserSet,'RealTime')
    UserSet.RealTime.FrametoRecon = UIValue;
else
    UserSet.FrametoRecon = UIValue;
end
GPUParams.NumFrame = UIValue;

%%% Check if 3D beamforming or 2D by checkking the name of the probe
if strcmpi(Trans.name(1),'M')
    sizeCoeff = [1,1,1];
else
    sizeCoeff = [1,1];
end

%%% Bmode/PI Realtime
switch UserSet.RealTime.Acquisition_Mode % Which kind of acquisition
    case 1 % BMode Transmission
        %%% Only one way to beamform, so we just have to create the
        %%% a matrix coefficient of ones having the shape of 
        %%% (numSample,numChannel,numAcq) 
        GPUParams.reconFramesCoefficient = ones([sizeCoeff GPUParams.numAcq]);

        numAcq_temp = length(GPUParams.reconFramesCoefficient);
        
    case 2 % PI Transmission
        %%% Two way to beamform:
        %%%     - Bmode (have to take care of pos/neg when we sum)
        %%%     - PI (we have less acquisition)
        switch UserSet.RealTime.Beamforming_Mode
        
            case 1 % Bmode Beamforming
                %%% Create a matrix coefficient of ones having the shape of 
                %%% (numSample,numChannel,numAcq) 
                GPUParams.reconFramesCoefficient = ones([sizeCoeff GPUParams.numAcq]);
                %%% we have a PI transmission but BMode beamforming we have to
                %%% provide a coefficient of -1 to every 2 acquisition othewise we will
                %%% sum NEG and POS acquisitions
                GPUParams.reconFramesCoefficient(1:2:end) = -1;
                numAcq_temp = length(GPUParams.reconFramesCoefficient);
                
            case 2 % PI Beamforming
                %%% we have a PI transmission and PI beamforming we have to
                %%% limit the beamforming to numAcq/2
                GPUParams.reconFramesCoefficient = ones([sizeCoeff GPUParams.numAcq/2]);
                numAcq_temp = length(GPUParams.reconFramesCoefficient);
        end        
end

% We use numAcq_temp to avoid to have too many frames for PI
GPUParams.reconFrames = ceil(linspace(1,numAcq_temp,GPUParams.NumFrame));
%%% As we don't beamforming everything, we just need the coefficient
%%% corresponding to the corresponding acquisition
GPUParams.reconFramesCoefficient = GPUParams.reconFramesCoefficient(GPUParams.reconFrames);

%%% Save the recon information in UserSet
    UserSet.RealTime.reconFrames = GPUParams.reconFrames;
    UserSet.RealTime.reconFramesCoefficient = GPUParams.reconFramesCoefficient;
    
%%% Need to update rf_data of GPUParams
GPUParams.rf_data=zeros(GPUParams.numSample,GPUParams.numChannel,...
    GPUParams.na,GPUParams.NumFrame,'single');
%%% Update initialization fo cuDAS_single
cuDAS_single(1,GPUParams.rf_data,GPUParams.degX_tx,GPUParams.degY_tx, ...
    GPUParams.actApertureX,GPUParams.actApertureY,GPUParams.actApertureZ,...
    GPUParams.filterCoef, GPUParams.pixelMapX,GPUParams.pixelMapY,GPUParams.pixelMapZ,...
    GPUParams.delay,GPUParams.fs, GPUParams.ftx,GPUParams.c, GPUParams.TXFocus,...
    GPUParams.SenseCutoff,GPUParams.ReconMode,GPUParams.gpuID);

assignin('base','UserSet',UserSet);
assignin('base','GPUParams',GPUParams);

end %CallBack-RTBeamformingFrameSlider

% - RTBeamformingProcessChoice
function RTBeamformingProcessChoice(~, ~, UIState)
% We choose which process (mean/std/max) we apply after the beamforming
% keyboard
UserSet = evalin('base','UserSet');
GPUParams = evalin('base','GPUParams');

if isfield(UserSet,'RealTime')
    UserSet.RealTime.ReconFrameProcess = UIState;
else
    UserSet.ReconFrameProcess = UIState;
end

GPUParams.ReconFrameProcess = UIState;

assignin('base','GPUParams',GPUParams);
assignin('base','UserSet',UserSet);
end %CallBack-RTBeamformingProcessChoice

% -SensCutoffSlider - Sensitivity cutoff change
function SensCutoffSlider(~, ~, UIValue)
GPUParams = evalin('base', 'GPUParams');
UserSet = evalin('base', 'UserSet');
GPUParams.SenseCutoff=UIValue;
 if isfield(UserSet,'RealTime')
     UserSet.RealTime.SenseCutoff = UIValue;
 else
     UserSet.SenseCutoff = UIValue;
 end
 
 if isfield(UserSet,'THI')
     UserSet.THI.SenseCutoff = UIValue;
 end
    
assignin('base','GPUParams',GPUParams);
assignin('base','UserSet',UserSet);

%%% UPDATE HERE, NEW BEAMFORMING OF CHEE HAU AND DIFFERENCE IF THI OR BMode
cuDAS_single(1,GPUParams.rf_data,GPUParams.degX_tx,GPUParams.degY_tx, ...
    GPUParams.actApertureX,GPUParams.actApertureY,GPUParams.actApertureZ,...
    GPUParams.filterCoef, GPUParams.pixelMapX,GPUParams.pixelMapY,GPUParams.pixelMapZ,...
    GPUParams.delay,GPUParams.fs, GPUParams.ftx,GPUParams.c, GPUParams.TXFocus,...
    GPUParams.SenseCutoff,GPUParams.ReconMode,GPUParams.gpuID);

return
end %CallBack-SensCutoffSlider

%CallBack-GPUDisplayCoefficientSlider
function GPUDisplayCoefficientSlider(~, ~, UIValue)
    GPUParams = evalin('base','GPUParams');
        GPUParams.GPURealTimeDisplayCoefficient = UIValue;
                
    UserSet=evalin('base','UserSet');

    %%% BMode or THI choice
    RealTimeDisplay_BModeTHI_choice = evalin('base','RealTimeDisplay_BModeTHI_choice');

    if RealTimeDisplay_BModeTHI_choice == 0 %0: B-Mode/PI
        if isfield(UserSet,'RealTime')
            UserSet.RealTime.GPURealTimeDisplayCoefficient = UIValue;
        else
            UserSet.GPURealTimeDisplayCoefficient = UIValue;
        end
    
    else %1: THI
        if isfield(UserSet,'THI')
            UserSet.THI.GPURealTimeDisplayCoefficient = UIValue;
        end
    end
    
    assignin('base','GPUParams',GPUParams);
    assignin('base','UserSet',UserSet);
end %CallBack-GPUDisplayCoefficientSlider

% -GPUFIlterCoefficientChoice - Choose GPU Filter coeffcieient
function GPUFIlterCoefficientChoice(~, ~,UIState)
% keyboard
GPUParams = evalin('base', 'GPUParams');
UserSet = evalin('base','UserSet');
Trans = evalin('base','Trans');

%%% In order to know which Filter we modify, BMode/PI and THI Mode may not
%%% have the same coefficient choice
RealTimeDisplay_BModeTHI_choice = evalin('base','RealTimeDisplay_BModeTHI_choice');
if RealTimeDisplay_BModeTHI_choice == 0 %1: Contrast Mode
    
    if isfield(UserSet,'RealTime')
        UserSet.RealTime.GPUFilterChoice = UIState;
        UserSetTemp = UserSet.RealTime;
    else
        UserSet.GPUFilterChoice = UIState;
        UserSetTemp = UserSet;
    end
        
else % THI Mode
    UserSet.THI.GPUFilterChoice = UIState;
    UserSetTemp = UserSet.THI;
end

if isfield(UserSetTemp,'GPUFilterCoeff')
    filterCoef = UserSetTemp.GPUFilterCoeff(UIState,:);
else
    filterCoef = Trans.Bandwidth;
end
disp(['New GPU Frequency Filter [' num2str(filterCoef(1)) ' ' num2str(filterCoef(2)) '] MHz']);

GPUParams.GPUFilterCoeff = filterCoef;
GPUParams.filterCoef= fir1(51,filterCoef/(Trans.frequency*2));
    clearvars filterCoef;

cuDAS_single(1,GPUParams.rf_data,GPUParams.degX_tx,GPUParams.degY_tx, ...
    GPUParams.actApertureX,GPUParams.actApertureY,GPUParams.actApertureZ,...
    GPUParams.filterCoef, GPUParams.pixelMapX,GPUParams.pixelMapY,GPUParams.pixelMapZ,...
    GPUParams.delay,GPUParams.fs, GPUParams.ftx,GPUParams.c, GPUParams.TXFocus,...
    GPUParams.SenseCutoff,GPUParams.ReconMode,GPUParams.gpuID);

    
assignin('base','GPUParams',GPUParams);
assignin('base','UserSet',UserSet);
return
end %CallBack-GPUFIlterCoefficientChoice

%CallBack-THIButton
function THIButton(~,~,~)
% keyboard
simMode = evalin('base','Resource.Parameters.simulateMode');
if simMode == 2
    return
end

UserSet = evalin('base','UserSet');

%%% Change colormap of the window, we need Resource
Resource = evalin('base','Resource');

%%% UIinterface
UIInterface = evalin('base','UIInterface');
%%% B-Mode/PI REAL-TIME event number
nTPCRealTime = evalin('base','nTPCRealTime');
%%% THI event number
nTPCTHI= evalin('base','nTPCTHI');
%%% BMode or THI choice
RealTimeDisplay_BModeTHI_choice = evalin('base','RealTimeDisplay_BModeTHI_choice');

%%% Slider of the voltage
hv1Sldr = findobj('Tag','hv1Sldr'); % Voltage currently used

    if RealTimeDisplay_BModeTHI_choice == 0 %0: B-Mode/PI to THI
        disp(['THI Imaging ' num2str(nTPCTHI)]);
        RealTimeDisplay_BModeTHI_choice = 1;
        nEvent = nTPCTHI;
%         %%% Change the Dynamic gain of the Process imageDisplay
%         evalin('base','Process(1).Parameters{8} = 5;'); % UPDATE 03/05/2019
            
        %%% Save the voltage currently used
            UserSet.RealTime.Voltage=get(hv1Sldr,'Value');
            
        %%% New voltage
            VOLTAGE = UserSet.THI.Voltage;
            
        %%% Take informaiton of the GPUParmas filter coeff
            GPUFilter_coeff = UserSet.THI.GPUFilterChoice;
        %%% Take information of the GPU Coefficient
            GPURealTimeDisplayCoefficient = UserSet.THI.GPURealTimeDisplayCoefficient;
        
        %%% We change the TGC to avoid to impact them
        TGC_choice = 2;
            
        Resource.DisplayWindow(1).Colormap = gray(256);
        
    else                      %1: THI mode to B-Mode/PI
        disp(['BMode/PI Imaging ' num2str(nTPCRealTime)]);
        RealTimeDisplay_BModeTHI_choice =0;
        nEvent = nTPCRealTime;
%         %%% Change the Dynamic gain of the Process imageDisplay
%         evalin('base','Process(1).Parameters{8} = 15;');
            
        %%% Save the voltage currently used
            UserSet.THI.Voltage=get(hv1Sldr,'Value');
            
        %%% New voltage
            VOLTAGE = UserSet.RealTime.Voltage;
            
        %%% Take informaiton of the GPUParmas filter coeff
            GPUFilter_coeff = UserSet.RealTime.GPUFilterChoice;
        %%% Take information of the GPU Coefficient
         GPURealTimeDisplayCoefficient = UserSet.RealTime.GPURealTimeDisplayCoefficient;
        
        %%% We change the TGC to avoid to impact them
            TGC_choice = 1;
            
        Resource.DisplayWindow(1).Colormap = pink(256);
    end
    assignin('base','RealTimeDisplay_BModeTHI_choice',RealTimeDisplay_BModeTHI_choice);
    assignin('base','UserSet',UserSet);       
         
     %%% UDPATE FOR COLORMAP
    newMap = Resource.DisplayWindow(1).Colormap;
    assignin('base','newMap',newMap);
    evalin('base','Resource.DisplayWindow(1).Colormap = newMap;');
    
    %%% Check if voltage VOLTAGE for TPC 1 (slider 1) is correct
    [result,hv] = setTpcProfileHighVoltage(VOLTAGE,1);
    %%% Change the voltage of the slider
        hv1Sldr = findobj('Tag','hv1Sldr'); % Voltage currently used
        set(hv1Sldr,'Value',hv);
    %%% Change the voltage text of the slider
        hv1Value = findobj('Tag','hv1Value'); % Voltage currently used
        set(hv1Value,'String',num2str(hv,'%.1f'));
        
    %%% Change the GPU Freq interface to show the right coeff
    %%% Frequency Contrast used
        index_interface = 6;
    RadioButton_GPUFilter = findobj(UIInterface{index_interface,1},UIInterface{index_interface,2});
    %%% Children value of the radiobutton set to 0
    RadioButton_Frequency_choice = zeros(1,length(RadioButton_GPUFilter.Children));
    %%% Set to 1 the GPUFilter coeff
    RadioButton_Frequency_choice(GPUFilter_coeff) = 1;
    %%% Radio button order is flipped in the handle
    RadioButton_Frequency_choice = fliplr(RadioButton_Frequency_choice);
    %%% Find which one is equal to 1
    RadioButton_Frequency_choice = find(RadioButton_Frequency_choice==1);
        %%% Apply the value to the parent handle, it changes automatically
        %%% the other one
        RadioButton_GPUFilter.Children(RadioButton_Frequency_choice).Value = 1;
                   
    %%% Change the GPURealTimeDisplayCoefficient interface to show the
    %%% right coeff
    %%% Handle of the slider
        index_interface = 5;
    Slider_GPURealTimeDisplayCoefficient = findobj(UIInterface{index_interface,1},UIInterface{index_interface,2});
    %%% Set Value of the Slider
    Slider_GPURealTimeDisplayCoefficient.Value = GPURealTimeDisplayCoefficient;
    %%% handle of the TextBox below the slider
    Slider_GPURealTimeDisplayCoefficient = findobj(UIInterface{index_interface,1},UIInterface{index_interface,4});
    %%% Set Value of the TextBox
    Slider_GPURealTimeDisplayCoefficient.String = num2str(GPURealTimeDisplayCoefficient); 
    
    %%% UPDATE FOR TPC
    %%% TGC handle
    test = findobj('Tag','TGCnum');
    %%% Find here how to control all the TGC and how to update all of them
    %%% in one shot
    set(test,'Value',TGC_choice);
    evalin('base',['if VDAS==1, Result = loadTgcWaveform(' num2str(TGC_choice) '); end']);
    
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Parameters',1,'startEvent',nEvent};
%%% Change colormap
    Control(2).Command = 'set&Run';
    Control(2).Parameters = {'DisplayWindow',1,'colormap',newMap};
%%% Change Voltage of the TPC 1 with the succeed value hv
    Control(3).Command = 'set&Run';
    Control(3).Parameters = {'TPC', 1, 'hv', hv};
    
%     Control(2).Command = 'set&Run';
%     Control(2).Parameters = {'TPC', TGC_choice};
evalin('base',['Resource.Parameters.startEvent =',num2str(nEvent),';'])

assignin('base','Control',Control);
end %CallBack-THIButton

%CallBack-AcquisitionButton - Jump to acquisition
function AcquisitionButton(~,~,~)
% keyboard
simMode = evalin('base','Resource.Parameters.simulateMode');
if simMode == 2
    return
end

%%% The Acquisition buttoin is updated to:
%%% - Generate a fodler where the acquisition will be saved
%%% - Modify the Process of Real Time saving to 

%%% Folder generated is based on the day
saving_buffer_counter = 1;
    assignin('base','saving_buffer_counter',saving_buffer_counter);
savedir = evalin('base','savedir');
        time_save = datestr(now,'HHMMSS');
        folder_acqui = [savedir 'CardiacGEM5ScD\' time_save '\'];
        mkdir(folder_acqui);
        assignin('base','folder_acqui',folder_acqui);
        assignin('base','time_save',time_save);
        
    %%% Temporal position for the name of the file
    UserSet=evalin('base','UserSet');
    filename_temp = ['PW_All_' num2str(UserSet.RealTime.NumFrame*UserSet.RealTime.numAcq) '_frames_'];
    %%% Filename generated is based on the day and time
    filename_acqui = [datestr(now,'yymmdd_HHMMSS') '_' filename_temp];
    assignin('base','filename_acqui',filename_acqui);
%     filename_acqui = 'TEST';
    
%%% Change the value of the Display to be sure we will initialize the right
%%% GPU beamforming
RealTimeDisplay_BModeTHI_choice = evalin('base','RealTimeDisplay_BModeTHI_choice');
RealTimeDisplay_BModeTHI_choice = 0;

assignin('base','RealTimeDisplay_BModeTHI_choice',RealTimeDisplay_BModeTHI_choice);
nAcquisition = evalin('base','nAcquisition');

Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',nAcquisition};
evalin('base',['Resource.Parameters.startEvent =',num2str(nAcquisition),';'])

assignin('base','Control',Control);
% assignin('base','freeze',1);
end %CallBack-AcquisitionButton

%% External functions with process object
%%% UPDATE Vantage V4.4

function RAZParameters() %Nom de la fonction
%-EF#1

    nTPCRealTime = evalin('base','nTPCRealTime');

    %%% For run with the right event
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Parameters',1,'startEvent',nTPCRealTime};
    evalin('base',['Resource.Parameters.startEvent =',num2str(nTPCRealTime),';'])
    assignin('base','Control',Control);
    % assignin('base','freeze',1);
    return
    
end %-EF#1

function saveRFfastData_Acquisition(RcvData) % Name of the function / Data comes from Receive Buffer 3 
%-EF#2

    % keyboard;
    
   
    saving_buffer_counter = evalin('base','saving_buffer_counter');
    
    % UserSet information
    UserSet=evalin('base','UserSet');

    % Folder where to save
    savedir = evalin('base','folder_acqui');
    time_save = evalin('base','time_save');
%     savedir = evalin('base','savedir');
%         time_save = datestr(now,'HHMMSS');
%         savedir = [savedir 'CardiacGEM5ScD\' time_save '\'];
%         mkdir(savedir);
    fname = [time_save,'_PW_ALL',['_' num2str(UserSet.RealTime.NumFrame*UserSet.RealTime.numAcq) '_frames']];
            assignin('base','savedir_THI',savedir);
            assignin('base','fname_THI',fname);
    GPUParams = evalin('base','GPUParams');
    tic;
    RcvData = RcvData(:,GPUParams.channelIndOrg,:); 
    RcvLastFrame=size(RcvData,3);

%     % Related to the size of the data, we may have to split it
%     S = whos('Data');
%     splitPart=ceil(S.bytes/2^32);
%      while rem(UserSet.RealTime.NumFrame/NumBuffer,splitPart)
%          splitPart = splitPart+1;
%      end
%      t = toc;
%      disp(['Start saving after ',num2str(t),'s']);
% 
%     tic
%     for indice = 1:splitPart
%         pos = (1:UserSet.RealTime.NumFrame/NumBuffer/splitPart)+(UserSet.RealTime.NumFrame/NumBuffer/splitPart * (indice-1));
%         RcvData = (Data(:,:,pos));
        savefast([savedir,fname '_Acquisition_' num2str(saving_buffer_counter)],'RcvData');
%     end
    t = toc;
    disp(['Buffer: ' num2str(saving_buffer_counter) ' - Acquisition saving in ',num2str(t),'s']);

    saving_buffer_counter = saving_buffer_counter+1;
    assignin('base','saving_buffer_counter',saving_buffer_counter);
    
    return
end %-EF#2

function GPU_init()
%-EF#3
% keyboard
%%% Last 2D/3D GPU Beamforming of Chee Hau Leow
GPUParams=struct('rf_data',[],'degX_tx',[],'degY_tx',[],...
    'actApertureX',[],'actApertureY',[],'actApertureZ',[],'filterCoef',[],...
    'pixelMapX',[],'pixelMapY',[], 'pixelMapZ',[],...
    'fs',[],'c',[],'delay',[],'TXFocus',[],...
    'SenseCutoff',0,'ReconMode',0,'gpuID',0,'channelInd',[]);

%Initialise GPUParams
Intermediate=evalin('base','Intermediate');
UserSet=evalin('base','UserSet');
RealTimeDisplay_BModeTHI_choice = evalin('base','RealTimeDisplay_BModeTHI_choice');

switch RealTimeDisplay_BModeTHI_choice
    
    case 0 % BMode
        %%% TOULEMONDE - Different definition of UserSet 
        if isfield(UserSet,'RealTime')
            UserSet = UserSet.RealTime;
        end
    case 1 % THI
%         keyboard
        %%% TOULEMONDE - Different definition of UserSet 
        if isfield(UserSet,'THI')
            UserSet = UserSet.THI;
        end
        %%% LRI Beamforming for THI
        GPUParams.ReconMode = 1;
end
% GPUParams.ReconMode = 1;

Resource=evalin('base','Resource');
Trans=evalin('base','Trans');
PData= evalin('base','PData');
%%%% TOULEMONDE - PData has two dimensions
    if size(PData,2)>1
       PData = PData(1);
    end
%%% TOULEMONDE - P
P=evalin('base','P');

%%% TOULEMONDE - Type of data saved (1) SUM - (2) PI - (3) AM
if isfield(UserSet,'TYPE')
    GPUParams.TYPE = UserSet.TYPE;
else
    GPUParams.TYPE = 1;
end

%%% TOULEMONDE - Type of aquisition and beamforming needed
if isfield(UserSet,'Acquisition_Mode')
    GPUParams.Acquisition_Mode = UserSet.Acquisition_Mode;
end
if isfield(UserSet,'Beamforming_Mode')
    GPUParams.Beamforming_Mode = UserSet.Beamforming_Mode;
end

GPUParams.SenseCutoff = UserSet.SenseCutoff;

%%% TOULEMONDE - Receive Apodization applied - Need to be int16
if isfield(UserSet,'RXApod')
    GPUParams.RXApod = UserSet.RXApod;
else
    GPUParams.RXApod = 1;
end

%%% UPDATE 13/01/2021 - Check if the ApertureTX exist
%%% TOULEMONDE - Number of Active aperture used in transmission
if isfield(UserSet,'ApertureTX') && isfield(UserSet.ApertureTX,'value')
    GPUParams.ApertureTX = UserSet.ApertureTX;
else
    GPUParams.ApertureTX.value = 1;
    GPUParams.ApertureTX.position = 1;
end

GPUParams.GPURealTimeDisplayCoefficient = UserSet.GPURealTimeDisplayCoefficient;

GPUParams.numSample=Intermediate.numRcvSamples;
GPUParams.numAcq=UserSet.numAcq;
GPUParams.na=UserSet.na;
%%% Update here about which frames to beamform
%%% We have defined previously which frames to beamform with which
%%% coefficient
if isfield(UserSet,'reconFramesCoefficient')
    GPUParams.NumFrame = UserSet.FrametoRecon;
    GPUParams.reconFramesCoefficient = UserSet.reconFramesCoefficient;
    GPUParams.reconFrames = UserSet.reconFrames;
else
    % Normal way of Chee Hau
    GPUParams.NumFrame=max(floor(UserSet.numAcq/UserSet.ReconSkipFrame),1);
    GPUParams.reconSkipFrame =UserSet.ReconSkipFrame;
    GPUParams.reconFramesCoefficient = 1;
    GPUParams.reconFrames= 1:GPUParams.reconSkipFrame:GPUParams.numAcq;
end

%%% Processing used during the beamforming when you beamform several frames
%%% ReconFrameProcess: 1:mean / 2:std / 3:max
if isfield(UserSet,'ReconFrameProcess')
    GPUParams.ReconFrameProcess = UserSet.ReconFrameProcess;
else
    GPUParams.ReconFrameProcess = 1; %default mean
end

if Trans.numelements<Resource.Parameters.numRcvChannels
    GPUParams.numChannel=Trans.numelements;
    GPUParams.channelInd = Trans.Connector;
else
    GPUParams.numChannel=Resource.Parameters.numRcvChannels;
    GPUParams.channelInd = Trans.Connector(1:GPUParams.numChannel);
end

%%% UPDATE 10/12/2021 GE Cardiac
%%% TOULEMONDE - Add name of Transducer
GPUParams.name = Trans.name;
%%% TOULEMONDE - Check GE Cardiac
 if strcmpi(GPUParams.name,'GEM5ScD')
    GPUParams.numChannel = 80;
    if isfield(UserSet,'GPUGEFull')
        GPUParams.GPUGEFull = UserSet.GPUGEFull;
    end
    GPUParams.channelIndOrg = GPUParams.channelInd; % save to have all channels
    GPUParams.channelInd = GPUParams.channelInd(1:80);
 else
    GPUParams.GPUGEFull = 0;
 end

%%% UPDATE Save own many channels received
GPUParams.colsPerFrame = Resource.RcvBuffer(1).colsPerFrame;

% GPUParams.rf_data=zeros(GPUParams.numSample,GPUParams.numChannel,GPUParams.na,GPUParams.NumFrame,'int16');
GPUParams.rf_data=zeros(GPUParams.numSample,GPUParams.numChannel,GPUParams.na,GPUParams.NumFrame,'single');

GPUParams.degX_tx(1)=0;
GPUParams.degY_tx(1)=0;
if UserSet.na>1
    for i=1:UserSet.na
        GPUParams.degX_tx(i)= UserSet.Angles(i);
%         (-UserSet.AngleRange/2+(i-1)*(UserSet.AngleRange)/(UserSet.na-1))*pi/180;
        GPUParams.degY_tx(i)=0;
    end
end

%%% Check if we have an acquisition with one or several apetures.
if GPUParams.ApertureTX.value > 1 || ((GPUParams.ApertureTX.value ==1) && (GPUParams.ApertureTX.position ~=1))   
    %%% Check resource, the number of channel used
    %%% If we have more than one aperture, we hve to see if there is an
    %%% overlap of the transmitted aperture 
    
    %%% Loop along the ApertureTX to obtain the connector position
    GPUParams.channelInd = zeros(GPUParams.numChannel,GPUParams.ApertureTX.value);
    for index_aperture = 1:GPUParams.ApertureTX.value
        GPUParams.channelInd(:,index_aperture) = Trans.Connector(Trans.HVMux.Aperture(:,GPUParams.ApertureTX.position(index_aperture))>0);
    end
    
    %%% If several apertures used, we need to know how many element are
    %%% overlapping and see how channels need to be used in reception for beamforming
    if GPUParams.ApertureTX.value > 1
        %%% Check how many elements are overlaping so how many needs to be
        %%% removed 
        
        %%% Initialyze the matrix we will mutiply the RcvData before
        %%% beamforming
        GPUParams.ApertureTX.RXApertureOverlap = ...
            zeros(GPUParams.numChannel,GPUParams.ApertureTX.value);
        
        Aperture = ...
            Trans.HVMux.ApertureES(:,GPUParams.ApertureTX.position);
        %%% Check between two consecutives aperture
        Aperture_Similar = logical(Aperture(:,1).*Aperture(:,2));
            [index] = find(Aperture_Similar==1);
        %%% - For Two apertures, the overlap is spread equivalently between
        %%% both aerptures
        %%% - For Three apertures, the central aperture have more element
        %%% Overlap Aperture 1/2: [1/4 3/4] overlap
        %%% Overlap Aperture 2/3: [3/4 1/4] overlap
        %%% Ex: if there are 64 elements overlaping between Aperture 1 and 2
        %%% 16 will be used for Aperture 1 and 48 for Aperture 2
        
        %%%% CHECK HERE
        switch GPUParams.ApertureTX.value
            case 2
                % For two apertures, we split the overlap between both
                % apertures.
                %%% Aperture 1
                    %%% Position of the first aperture
                    pos = [1:GPUParams.numChannel-round(length(index)/2)];
                %%% Set all element of aperture equal to 1
                GPUParams.ApertureTX.RXApertureOverlap(pos,1) = 1;
                %%% Find the index of the aperture where it is equal to
                %%% one in the 128 aperture
                GPUParams.ApertureTX.RXApertureOverlap_index{:,1} = find(GPUParams.ApertureTX.RXApertureOverlap(:,1)==1);
                %%% Find the first element of the aperture to when we
                %%% will copy the 128 to the 256 aperture to set them
                %%% at the right position
                GPUParams.ApertureTX.RXApertureOverlap_index_min(:,1) = single(min(GPUParams.ApertureTX.RXApertureOverlap_index{:,1})-1);
                %%% The element in the aperture are not from 1-128 but
                %%% are corresponding to the pin of the VErasonics and
                %%% it is probe dependent. We check in channelInd to
                %%% use the right pins for the aperture RXApertureOverlap_index
                GPUParams.ApertureTX.RXApertureOverlap_index{:,1} = GPUParams.channelInd(GPUParams.ApertureTX.RXApertureOverlap_index{:,1},1);
                
                %%% Aperture 2
                    %%% Position of the second aperture
                    pos = [round(length(index)/2):GPUParams.numChannel];
                %%% Set all element of aperture equal to 1
                GPUParams.ApertureTX.RXApertureOverlap(pos,2) = 1;
                %%% Find the index of the aperture where it is equal to
                %%% one in the 128 aperture
                GPUParams.ApertureTX.RXApertureOverlap_index{:,2} = find(GPUParams.ApertureTX.RXApertureOverlap(:,2)==1);
                %%% Find the first element of the aperture to when we
                %%% will copy the 128 to the 256 aperture to set them
                %%% at the right position
                GPUParams.ApertureTX.RXApertureOverlap_index_min(:,2) = single(min(GPUParams.ApertureTX.RXApertureOverlap_index{:,2})-1);
                %%% The element in the aperture are not from 1-128 but
                %%% are corresponding to the pin of the VErasonics and
                %%% it is probe dependent. We check in channelInd to
                %%% use the right pins for the aperture RXApertureOverlap_index
                GPUParams.ApertureTX.RXApertureOverlap_index{:,2} = GPUParams.channelInd(GPUParams.ApertureTX.RXApertureOverlap_index{:,2},2);
                
            case 3
                % For three apertures, the overlap is different depending
                % of the aperture
                %%% Aperture 1
                    %%% Position of the first aperture
                    pos = [1:GPUParams.numChannel-round(length(index)*3.0/4)];
                %%% Set all element of aperture equal to 1
                GPUParams.ApertureTX.RXApertureOverlap(pos,1) = 1;
                %%% Find the index of the aperture where it is equal to
                %%% one in the 128 aperture
                GPUParams.ApertureTX.RXApertureOverlap_index{:,1} = find(GPUParams.ApertureTX.RXApertureOverlap(:,1)==1);
                %%% Find the first element of the aperture to when we
                %%% will copy the 128 to the 256 aperture to set them
                %%% at the right position
                GPUParams.ApertureTX.RXApertureOverlap_index_min(:,1) = single(min(GPUParams.ApertureTX.RXApertureOverlap_index{:,1})-1);
                %%% The element in the aperture are not from 1-128 but
                %%% are corresponding to the pin of the VErasonics and
                %%% it is probe dependent. We check in channelInd to
                %%% use the right pins for the aperture RXApertureOverlap_index
                GPUParams.ApertureTX.RXApertureOverlap_index{:,1} = GPUParams.channelInd(GPUParams.ApertureTX.RXApertureOverlap_index{:,1},1);
                
                %%% Aperture 2
                    %%% Position of the second aperture
                    pos = [round(length(index)*1.0/4):GPUParams.numChannel-round(length(index)*1.0/4)-1];
                %%% Set all element of aperture equal to 1
                GPUParams.ApertureTX.RXApertureOverlap(pos,2) = 1;
                %%% Find the index of the aperture where it is equal to
                %%% one in the 128 aperture
                GPUParams.ApertureTX.RXApertureOverlap_index{:,2} = find(GPUParams.ApertureTX.RXApertureOverlap(:,2)==1);
                %%% Find the first element of the aperture to when we
                %%% will copy the 128 to the 256 aperture to set them
                %%% at the right position
                GPUParams.ApertureTX.RXApertureOverlap_index_min(:,2) = single(min(GPUParams.ApertureTX.RXApertureOverlap_index{:,2})-1);
                %%% The element in the aperture are not from 1-128 but
                %%% are corresponding to the pin of the VErasonics and
                %%% it is probe dependent. We check in channelInd to
                %%% use the right pins for the aperture RXApertureOverlap_index
                GPUParams.ApertureTX.RXApertureOverlap_index{:,2} = GPUParams.channelInd(GPUParams.ApertureTX.RXApertureOverlap_index{:,2},2);
                
                %%% Aperture 3
                    %%% Position of the third aperture
                    pos = [round(length(index)*3.0/4):GPUParams.numChannel];
                %%% Set all element of aperture equal to 1
                GPUParams.ApertureTX.RXApertureOverlap(pos,3) = 1;
                %%% Find the index of the aperture where it is equal to
                %%% one in the 128 aperture
                GPUParams.ApertureTX.RXApertureOverlap_index{:,3} = find(GPUParams.ApertureTX.RXApertureOverlap(:,3)==1);
                %%% Find the first element of the aperture to when we
                %%% will copy the 128 to the 256 aperture to set them
                %%% at the right position
                GPUParams.ApertureTX.RXApertureOverlap_index_min(:,3) = single(min(GPUParams.ApertureTX.RXApertureOverlap_index{:,3})-1);
                %%% The element in the aperture are not from 1-128 but
                %%% are corresponding to the pin of the VErasonics and
                %%% it is probe dependent. We check in channelInd to
                %%% use the right pins for the aperture RXApertureOverlap_index
                GPUParams.ApertureTX.RXApertureOverlap_index{:,3} = GPUParams.channelInd(GPUParams.ApertureTX.RXApertureOverlap_index{:,3},3);
        end
        clearvars Aperture;
        clearvars Aperture_Similar index pos;
    else
        %%% We set the aperture overlap to 1
        GPUParams.ApertureTX.RXApertureOverlap = ones(GPUParams.numChannel,1);        
        GPUParams.ApertureTX.RXApertureOverlap_index_min(:,1) = 1;
        GPUParams.ApertureTX.RXApertureOverlap_index{:,1} = find(GPUParams.ApertureTX.RXApertureOverlap(:,1)==1);
        GPUParams.ApertureTX.RXApertureOverlap_index{:,1} = GPUParams.channelInd(GPUParams.ApertureTX.RXApertureOverlap_index{:,1},1);
    end
    %%% Size of the beamforming will be bigger as we want to use all
    %%% elements
    GPUParams.numChannel = length(Trans.Connector);
    GPUParams.rf_data=zeros(GPUParams.numSample,GPUParams.numChannel,GPUParams.na,GPUParams.NumFrame,'single');
else
    %%% Even if we have no aperture used, we initialize RXApertureOverlap to
    %%% multiply our RcvData before beamforming
    GPUParams.ApertureTX.RXApertureOverlap = ones(GPUParams.numChannel,1);        
    GPUParams.ApertureTX.RXApertureOverlap_index_min(:,1) = 1;
    GPUParams.ApertureTX.RXApertureOverlap_index{:,1} = find(GPUParams.ApertureTX.RXApertureOverlap(:,1)==1);
    GPUParams.ApertureTX.RXApertureOverlap_index{:,1} = GPUParams.channelInd(GPUParams.ApertureTX.RXApertureOverlap_index{:,1},1);
end

%%% Define active aperture of the probe
GPUParams.actApertureX=Trans.ElementPos(:,1)*1e-3;
GPUParams.actApertureY=Trans.ElementPos(:,2)*1e-3;
GPUParams.actApertureZ=Trans.ElementPos(:,3)*1e-3;

%%% TOULEMONDE - Provide filter coefficient for GPU Beamforming
if isfield(UserSet,'GPUFilterCoeff')
    filterCoef = UserSet.GPUFilterCoeff(UserSet.GPUFilterChoice,:);
else
    filterCoef = Trans.Bandwidth;
end
GPUParams.GPUFilterCoeff = filterCoef;
GPUParams.filterCoef= fir1(51,filterCoef/(Trans.frequency*2));
    clearvars filterCoef;
    
GPUParams.fs=Trans.frequency*1e6*4;
GPUParams.ftx =UserSet.TXFreq*1e6;
GPUParams.c = Resource.Parameters.speedOfSound;
GPUParams.TXFocus =UserSet.TXFocusMM*1e-3;
GPUParams.delay=(2*Trans.lensCorrection*1e-3/GPUParams.c)*ones(UserSet.na,1);

% GPUParams.pixelMapX = linspace(min(Trans.ElementPos(:)),max(Trans.ElementPos(:)),PData.Size(2))*1e-3;
GPUParams.pixelMapX = -linspace(PData.Origin(1),-PData.Origin(1),PData.Size(2))*GPUParams.c/(Trans.frequency*1e6); %-(UserSet.aperture*Trans.spacingMm/1e3)
GPUParams.pixelMapY=0;

%%% TRIAL RSVD - UPDATE
if ~isfield(UserSet,'SVDEn_Random')
    GPUParams.SVDEn_Random = 50;
else
    GPUParams.SVDEn_Random = UserSet.SVDEn_Random;
end
    svd_img_data = single(zeros(PData.Size(1),PData.Size(2),GPUParams.SVDEn_Random));

wavelength=GPUParams.c/Trans.frequency*1e-6;
% range = evalin('base','UIParam.Range');
%%% TOULEMONDE UDAPTE, range is realated to P.endDpeth
range = P(1).endDepth;
GPUParams.pixelMapZ = linspace(0,range*wavelength,PData.Size(1));
% GPUParams.pixelMapZ = linspace(0,UserSet.sampleDepthMM,PData.Size(1))*1e-3;

GPUParams.IQData= zeros(PData.Size(1),PData.Size(2),GPUParams.NumFrame,'single'); % ICI
GPUParams.IQData= complex(GPUParams.IQData,0);
GPUParams.SVDEn = UserSet.SVDEn; 

cuDAS_single(1,GPUParams.rf_data,GPUParams.degX_tx,GPUParams.degY_tx, ...
    GPUParams.actApertureX,GPUParams.actApertureY,GPUParams.actApertureZ,...
    GPUParams.filterCoef, GPUParams.pixelMapX,GPUParams.pixelMapY,GPUParams.pixelMapZ,...
    GPUParams.delay,GPUParams.fs, GPUParams.ftx,GPUParams.c, GPUParams.TXFocus,...
    GPUParams.SenseCutoff,GPUParams.ReconMode,GPUParams.gpuID);
        
assignin('base','GPUParams',GPUParams);
assignin('base','svd_img_data',svd_img_data);

return
end %-EF#3

function img =GPU_beamforming(RcvData)
%-EF#4
% keyboard;
tic;

GPUParams=evalin('base','GPUParams');
% sizeRF = size(GPUParams.rf_data);

%%% ORG
% % temp=permute(reshape(RcvData(1:(GPUParams.numSample*GPUParams.na*GPUParams.numAcq),:,:)...
% %     ,GPUParams.numSample,GPUParams.na,GPUParams.numAcq,[]),[1,4,2,3]);

%%% Check if it is BMode or PI acquisition 
%%% If it is PI or BMode beamforming
%%% How we process them (mean / std / max)
%%% See with Kai how we can accumulate the previous beamform (an option
%%% from Verasonics certainly)

%%% NEW
%%% Update 13/01/2021 - Take into account the number of aperture
Start_DATA = 0;
Length_DATA = (GPUParams.numSample*GPUParams.na*GPUParams.numAcq*GPUParams.ApertureTX.value);
    
%%% Update 10/12/2021
%%% Check if it is the GEM5ScD in order to use only 160 or 80 channels
 if strcmpi(GPUParams.name,'GEM5ScD')
     RcvData = RcvData(:,GPUParams.channelIndOrg,:);
     if GPUParams.GPUGEFull
         RcvData = RcvData(:,1:80,:) + RcvData(:,81:160,:);
     else
         RcvData = RcvData(:,1:80,:);
     end
     GPUParams.channelInd = 1:80;
 end
 
if isfield(GPUParams,'Beamforming_Mode')
    %%% PI Beamforming
    if GPUParams.Beamforming_Mode == 2
        %%% Devide by two the length of the acquisition
        Length_DATA = Length_DATA/2;
        %%% Divide by two the numAcq here as it will not be exported in
        %%% order to do a good unwrap
        GPUParams.numAcq = GPUParams.numAcq/2;        
        %%% Change he TYPe into 4 as PI to do the SUM
        %%% We have create a new GPU Case as NEG and POS are next each
        %%% other and not all NEG then all POS in the memory:
        %%% Normal case in memory for 1 frame [A1N/A2N/A3N....A1P/A2P/A3P]
        %%%  This case in memory for 1 frame  [A1N/A1P/A2N/A2P/A3N/A3P...]
        GPUParams.TYPE = 2;
    end
end

%%% We have to adapt depending of the arealdy SUMMED, PI or AM
switch GPUParams.TYPE
    case 1 % Memory already SUM
        POSITION1 = Start_DATA+[1:Length_DATA];
        %%% Update 13/01/2021 - The order of the element will be done later
%         temp = RcvData(POSITION1,GPUParams.channelInd,:);
        temp = RcvData(POSITION1,:,:);
        
    case 2 % We have Positive and Negative Data not summed so we have to extract and sum them
        
         %%% Normal case in memory for 1 frame [A1N/A2N/A3N....A1P/A2P/A3P]
        POSITION1 = Start_DATA+[1:Length_DATA];
        POSITION2 = (Start_DATA+Length_DATA)+[1:Length_DATA];
        
        %%%  This case in memory for 1 frame  [A1N/A1P/A2N/A2P/A3N/A3P...]
        % We have Negative and POsitive not summed and they are following
        if isfield(GPUParams,'Beamforming_Mode')
            %%% PI Beamforming
            if GPUParams.Beamforming_Mode == 2
                Length_DATA = GPUParams.numSample;
                index_neg = [0:2:(GPUParams.numAcq-1)*2];
                
                POSITION1 = repmat([1:Length_DATA]',[1 GPUParams.numAcq]);
                POSITION1 = POSITION1+index_neg*Length_DATA;
                
                POSITION2 = POSITION1+Length_DATA;
            end
        end
        
        %%% Update 13/01/2021 - The order of the element will be done later
%         temp = RcvData(POSITION1,GPUParams.channelInd,:)+RcvData(POSITION2,GPUParams.channelInd,:);
        temp = RcvData(POSITION1,:,:)+RcvData(POSITION2,:,:);
        
    case 3 % We have Half / Full / Half Data not summed so we have to extract and sum them
        POSITION1 = Start_DATA+[1:Length_DATA];
        POSITION2 = (Start_DATA+Length_DATA)+[1:Length_DATA];
        POSITION3 = (Start_DATA+2*Length_DATA)+[1:Length_DATA];
        %%% Update 13/01/2021 - The order of the element will be done later
%         temp = RcvData(POSITION1,GPUParams.channelInd,:)+RcvData(POSITION2,GPUParams.channelInd,:)+RcvData(POSITION3,GPUParams.channelInd,:);
        temp = RcvData(POSITION1,:,:) + RcvData(POSITION2,:,:) + RcvData(POSITION3,:,:);

end

%%% Reshape to have the data: [depth,angle,aperture,numAcq,channel]
temp=reshape(temp,GPUParams.numSample,GPUParams.na,GPUParams.ApertureTX.value,GPUParams.numAcq,[length(GPUParams.channelInd)]);
%%% Permute to have the data: [depth,channel,angle,numAcq,aperture]
temp = permute(temp,[1,5,2,4,3]);

%%% RXApod is moved
% temp = temp.*(GPUParams.RXApod);
   
for index_aperture=1:GPUParams.ApertureTX.value
    if GPUParams.ApertureTX.value > 1 || ((GPUParams.ApertureTX.value ==1) && (GPUParams.ApertureTX.position ~=1))
%         GPUParams.rf_data = GPUParams.rf_data*0;
    %%% If we have a multiplexed probe, we create GPUParams.aperture    
%         GPUParams.rf_data(:,GPUParams.ApertureTX.position(index_aperture)+[0:length(GPUParams.channelInd)-1],:,:) = temp(:,:,:,GPUParams.reconFrames,index_aperture).*GPUParams.RXApod;
        GPUParams.rf_data(:,GPUParams.ApertureTX.position(index_aperture)+[0:length(GPUParams.ApertureTX.RXApertureOverlap_index{:,index_aperture})-1]+GPUParams.ApertureTX.RXApertureOverlap_index_min(:,index_aperture),:,:) = ...
            single(temp(:,GPUParams.ApertureTX.RXApertureOverlap_index{:,index_aperture},:,GPUParams.reconFrames,index_aperture));
    else
       %%% Convert to single and multiplication to apodization
       GPUParams.rf_data=single(temp(:,GPUParams.channelInd,:,GPUParams.reconFrames,index_aperture)).*GPUParams.RXApod; 
    end
end
%%%  Even if rf_data is not the same size as define in gpu_init we will
%%%  provide the GPUParams.reconFrames to have thr right size
% GPUParams.rf_data=permute(reshape(temp...
%     ,GPUParams.numSample,GPUParams.na,GPUParams.numAcq,[length(GPUParams.channelInd)]),[1,4,2,3]);

t_1 = toc;

%%% BEAMFORMING
% img = cuDAS_single(0,GPUParams.rf_data);
% GPUParams.rf_data(:,36,:) = 0;
img = cuDAS_single(0,GPUParams.rf_data);

% Multuply the img beamformed image by reconFramesCoefficient to apply
% coefficient (default 1)
GPUParams.IQData=img.*GPUParams.reconFramesCoefficient;

t_2 = toc-(t_1);

% disp(['Time - t1 : ' num2str(t_1,'%1.4f') ' t2 : ' num2str(t_2,'%1.4f')]);
    
%     keyboard
    %%% TOULEMONDE - IF WE WANT to AVERAGE / STD / MAX consecutives frames 
    switch GPUParams.ReconFrameProcess
            
        case 2 % STD
%             if GPUParams.ReconMode == 1
%                 b = nchoosek([1:GPUParams.na],2);
%                 test = img(:,:,b(:,1)).*img(:,:,b(:,2));
%                 test_sign = sign(test);
%                 img = abs(double(GPUParams.GPURealTimeDisplayCoefficient*sum(sqrt(abs(test)).*test_sign,3)));
%             else
                img=double(GPUParams.GPURealTimeDisplayCoefficient*std(abs(GPUParams.IQData),0,3));
%             end
        case 3 % SVD

            if GPUParams.SVDEn==1
                if size(GPUParams.IQData,3)==1
                   disp('SVD type 1 need at least 2 frames');
                   img=double(GPUParams.GPURealTimeDisplayCoefficient*mean(abs(GPUParams.IQData),3));
                   return
                end
                
                % Create Handle For Figure To Determine Cutoff
                persistent my_handleSVD

                if isempty(my_handleSVD)
                    figure('Name','SVD Chee Hau');
                    my_handleSVD = axes();
                    xlabel(my_handleSVD,'Singular Value');ylabel(my_handleSVD,'Energy [dB]');
                        title('Energy of Singular value');
                end
                set(my_handleSVD,'NextPlot','replacechildren');
                
                %Clutter Filtering
                sizeIQ=size(GPUParams.IQData);
                ns=sizeIQ(1)*sizeIQ(2);
                permInd =randperm(ns);

                for i=10:-1:1
                    if mod(ns,i)==0
                        nsplit = i;
                        break;
                    end
                end

                data=reshape(GPUParams.IQData,[],sizeIQ(3));
                permInd =reshape(permInd,[],nsplit);

                %%
                %%% HAVE TO UPDATE HERE WITH SVD AS A CHOICE AS 
                dataFilter = zeros(ns,sizeIQ(3),'like',data);
                for i=1:nsplit
            %         tic;
            %         disp(['SVD processing Part : ',num2str(i),'\',num2str(nsplit)])
                    tempData=data(permInd(:,i),:);
                    [U,S,V]=svd(tempData,'econ');
                    
            %         %Automatic choosing threshold
            %         logEnergy=20*log10(diag(S)/max(diag(S)));
            %         logEnergy1=sgolayfilt(double(logEnergy),2,21);
            %         %     logEnergy1=sgolayfilt(double(logEnergy),2,5);
            %         dy1=gradient(logEnergy1);
            %         dy2=gradient(dy1);
            %         curve=abs((dy2)./(1+dy1.^2).^(3/2));
            %         [~,index1]=max(curve(1:ceil(0.5*length(curve))));
            %         [~,index2]=max(curve(ceil(0.5*length(curve))+1:end));
            %         index2=index2+ceil(0.5*length(curve));

            % %         Test for distance
                    logEnergy=20*log10(diag(S)/max(diag(S)));
                    logEnergy1=logEnergy-min(logEnergy(:));
                    logEnergy1=logEnergy1/max(logEnergy1);
                    xlog=((1:length(logEnergy1))/length(logEnergy1))';
                    dist=sqrt(logEnergy1.^2+xlog.^2);
                    %              plot(xlog,logEnergy);
                    [~,index1]=min(dist);
                    index2 = sizeIQ(3);
            %         index1=max(index1,index1a);
                    dataFilter(permInd(:,i),:)=U*S(:,index1:end)*V(:,index1:end)';
            %         toc
            
                    cla(my_handleSVD);
                    hold(my_handleSVD,'all');
                    plot(my_handleSVD,logEnergy);
                    plot(my_handleSVD,index1,logEnergy(index1),'x');
                    plot(my_handleSVD,index2,logEnergy(index2),'x');
                    drawnow limitrate;
    
                end

                img=double(reshape(mean(dataFilter(:,1:end-1).*conj(dataFilter(:,2:end)),2),sizeIQ(1),sizeIQ(2)));
%                 img=double(reshape(dataFilter(:,:),sizeIQ(1),sizeIQ(2),sizeIQ(3)));
                img = sqrt(img);
                img = GPUParams.GPURealTimeDisplayCoefficient*mean(abs(img),3);
            %     img=double(GPUParams.GPURealTimeDisplayCoefficient*mean(abs(GPUParams.IQData),3));

            else

                % Create Handle For Figure To Determine Cutoff
                persistent my_handle

                if isempty(my_handle)
                    figure();
                    my_handle = axes();
                    set(my_handle,'Ydir','reverse');
                    axis(my_handle,'image');
                end
                set(my_handle,'NextPlot','replacechildren');

                % Load Data Stack
                svd_img_data = evalin('base','svd_img_data');

                % Load Instantaneous Data
                RData=GPUParams.IQData;

                % Image Size
                nx = size(RData,2);
                nz = size(RData,1);
                nt = GPUParams.SVDEn_Random;

                % Replace last frames with new frames
                svd_img_data(:,:,1:size(RData,3))  = []; 
                svd_img_data(:,:,nt-size(RData,3)+1:nt) = RData; 

                % SVD Filter RSVD
                k=2; % Singular Values you want to remove
                
                CasoratiData = reshape(svd_img_data,[nz*nx,nt]);
                [U,S,V] = SubFunctionRSVD(single(CasoratiData),k);
                SVDData = U*S*V';
                SVDData = reshape(SVDData,[nz nx nt]);

                % Plot Spatial Similarity between Singular Values
                if ~isempty(U) 
                    U = reshape(U,[nz nx]);
                    imagesc(corrcoef(abs(U)),'Parent',my_handle);
                    drawnow limitrate
                end

                % Save Image for Display
                img=double(GPUParams.GPURealTimeDisplayCoefficient*(abs(SVDData(:,:,end))));

                % Export Data to Verasonics
                assignin('base','svd_img_data',svd_img_data);
            end
            
        otherwise % MEAN by default
%             if exist('my_handleSVD','var')
            if GPUParams.ReconMode == 1
                img=double(GPUParams.GPURealTimeDisplayCoefficient*mean(abs(squeeze(sum(GPUParams.IQData,3))),3));
            else
                img=double(GPUParams.GPURealTimeDisplayCoefficient*mean(abs(GPUParams.IQData),3));
            end
    end

return
end %-EF#4

function img =GPU_beamforming_THI(RcvData)
%-EF#5
% Exactly the same script as previous but with another buffer
% keyboard;
tic;

GPUParams=evalin('base','GPUParams');
% sizeRF = size(GPUParams.rf_data);

%%% ORG
% % temp=permute(reshape(RcvData(1:(GPUParams.numSample*GPUParams.na*GPUParams.numAcq),:,:)...
% %     ,GPUParams.numSample,GPUParams.na,GPUParams.numAcq,[]),[1,4,2,3]);

%%% Check if it is BMode or PI acquisition 
%%% If it is PI or BMode beamforming
%%% How we process them (mean / std / max)
%%% See with Kai how we can accumulate the previous beamform (an option
%%% from Verasonics certainly)

%%% NEW
%%% Update 13/01/2021 - Take into account the number of aperture
Start_DATA = 0;
Length_DATA = (GPUParams.numSample*GPUParams.na*GPUParams.numAcq*GPUParams.ApertureTX.value);
   
%%% Update 10/12/2021
%%% Check if it is the GEM5ScD in order to use only 160 or 80 channels
 if strcmpi(GPUParams.name,'GEM5ScD')
     RcvData = RcvData(:,GPUParams.channelIndOrg,:);
     if GPUParams.GPUGEFull
         RcvData = RcvData(:,1:80,:) + RcvData(:,81:160,:);
     else
         RcvData = RcvData(:,1:80,:);
     end
     GPUParams.channelInd = 1:80;
 end
 
if isfield(GPUParams,'Beamforming_Mode')
    %%% PI Beamforming
    if GPUParams.Beamforming_Mode == 2
        %%% Devide by two the length of the acquisition
        Length_DATA = Length_DATA/2;
        %%% Divide by two the numAcq here as it will not be exported in
        %%% order to do a good unwrap
        GPUParams.numAcq = GPUParams.numAcq/2;        
        %%% Change he TYPe into 4 as PI to do the SUM
        %%% We have create a new GPU Case as NEG and POS are next each
        %%% other and not all NEG then all POS in the memory:
        %%% Normal case in memory for 1 frame [A1N/A2N/A3N....A1P/A2P/A3P]
        %%%  This case in memory for 1 frame  [A1N/A1P/A2N/A2P/A3N/A3P...]
        GPUParams.TYPE = 2;
    end
end

%%% We have to adapt depending of the arealdy SUMMED, PI or AM
switch GPUParams.TYPE
    case 1 % Memory already SUM
        POSITION1 = Start_DATA+[1:Length_DATA];
        %%% Update 13/01/2021 - The order of the element will be done later
%         temp = RcvData(POSITION1,GPUParams.channelInd,:);
        temp = RcvData(POSITION1,:,:);
        
    case 2 % We have Positive and Negative Data not summed so we have to extract and sum them
        
         %%% Normal case in memory for 1 frame [A1N/A2N/A3N....A1P/A2P/A3P]
        POSITION1 = Start_DATA+[1:Length_DATA];
        POSITION2 = (Start_DATA+Length_DATA)+[1:Length_DATA];
        
        %%%  This case in memory for 1 frame  [A1N/A1P/A2N/A2P/A3N/A3P...]
        % We have Negative and POsitive not summed and they are following
        if isfield(GPUParams,'Beamforming_Mode')
            %%% PI Beamforming
            if GPUParams.Beamforming_Mode == 2
                Length_DATA = GPUParams.numSample;
                index_neg = [0:2:(GPUParams.numAcq-1)*2];
                
                POSITION1 = repmat([1:Length_DATA]',[1 GPUParams.numAcq]);
                POSITION1 = POSITION1+index_neg*Length_DATA;
                
                POSITION2 = POSITION1+Length_DATA;
            end
        end
        
        %%% Update 13/01/2021 - The order of the element will be done later
%         temp = RcvData(POSITION1,GPUParams.channelInd,:)+RcvData(POSITION2,GPUParams.channelInd,:);
        temp = RcvData(POSITION1,:,:)+RcvData(POSITION2,:,:);
        
    case 3 % We have Half / Full / Half Data not summed so we have to extract and sum them
        POSITION1 = Start_DATA+[1:Length_DATA];
        POSITION2 = (Start_DATA+Length_DATA)+[1:Length_DATA];
        POSITION3 = (Start_DATA+2*Length_DATA)+[1:Length_DATA];
        %%% Update 13/01/2021 - The order of the element will be done later
%         temp = RcvData(POSITION1,GPUParams.channelInd,:)+RcvData(POSITION2,GPUParams.channelInd,:)+RcvData(POSITION3,GPUParams.channelInd,:);
        temp = RcvData(POSITION1,:,:) + RcvData(POSITION2,:,:) + RcvData(POSITION3,:,:);

end

%%% Reshape to have the data: [depth,angle,aperture,numAcq,channel]
temp=reshape(temp,GPUParams.numSample,GPUParams.na,GPUParams.ApertureTX.value,GPUParams.numAcq,[length(GPUParams.channelInd)]);
%%% Permute to have the data: [depth,channel,angle,numAcq,aperture]
temp = permute(temp,[1,5,2,4,3]);

for index_aperture=1:GPUParams.ApertureTX.value
    if GPUParams.ApertureTX.value > 1 || ((GPUParams.ApertureTX.value ==1) && (GPUParams.ApertureTX.position ~=1))
%         GPUParams.rf_data = GPUParams.rf_data*0;
    %%% If we hve a multiplexied probe, we create GPUParams.aperture    
%         GPUParams.rf_data(:,GPUParams.ApertureTX.position(index_aperture)+[0:length(GPUParams.channelInd)-1],:,:) = temp(:,:,:,GPUParams.reconFrames,index_aperture).*GPUParams.RXApod;
        GPUParams.rf_data(:,GPUParams.ApertureTX.position(index_aperture)+[0:length(GPUParams.ApertureTX.RXApertureOverlap_index{:,index_aperture})-1]+GPUParams.ApertureTX.RXApertureOverlap_index_min(:,index_aperture),:,:) = ...
            single(temp(:,GPUParams.ApertureTX.RXApertureOverlap_index{:,index_aperture},:,GPUParams.reconFrames,index_aperture));
    else
       %%% Convert to single and multiplication to apodization
       GPUParams.rf_data=single(temp(:,GPUParams.channelInd,:,GPUParams.reconFrames,index_aperture)).*GPUParams.RXApod; 
    end
end

% %%% If we hve a multiplexied probe, we create GPUParams.aperture
% if isfield(GPUParams,'aperture')    
%     GPUParams.rf_data(:,GPUParams.aperture+[0:length(GPUParams.channelInd)-1],:,:) = temp(:,:,:,GPUParams.reconFrames);
% else
%    GPUParams.rf_data=temp(:,:,:,GPUParams.reconFrames); 
% end

%%%  Even if rf_data is not the same size as define in gpu_init we will
%%%  provide the GPUParams.reconFrames to have thr right size
% GPUParams.rf_data=permute(reshape(temp...
%     ,GPUParams.numSample,GPUParams.na,GPUParams.numAcq,[length(GPUParams.channelInd)]),[1,4,2,3]);

t_1 = toc;

%%% BEAMFORMING


% GPUParams.rf_data(:,36,:) = 0;

img = cuDAS_single(0,GPUParams.rf_data);
% Multuply the img beamformed image by reconFramesCoefficient to apply
% coefficient (default 1)
if GPUParams.ReconMode == 1
    GPUParams.IQData=sum(img,3).*GPUParams.reconFramesCoefficient;
else
    GPUParams.IQData=img.*GPUParams.reconFramesCoefficient;
end

t_2 = toc-(t_1);

% disp(['Time - t1 : ' num2str(t_1,'%1.4f') ' t2 : ' num2str(t_2,'%1.4f')]);

    %%% TOULEMONDE - IF WE WANT to AVERAGE / STD / MAX consecutives frames 
    switch GPUParams.ReconFrameProcess
            
        case 2 % STD
            if GPUParams.ReconMode == 1
                b = nchoosek([1:GPUParams.na],2);
                test = img(:,:,b(:,1)).*img(:,:,b(:,2));
                test_sign = sign(test);
                img = abs(double(GPUParams.GPURealTimeDisplayCoefficient*sum(sqrt(abs(test)).*test_sign,3)));
%                 img = double(GPUParams.GPURealTimeDisplayCoefficient*sum(abs(img),3));
            else
                img=double(GPUParams.GPURealTimeDisplayCoefficient*std(abs(GPUParams.IQData),0,3));
            end
        case 3 % MAX
            img=double(GPUParams.GPURealTimeDisplayCoefficient*max(abs(GPUParams.IQData),[],3));
        otherwise % MEAN by default
            img=double(GPUParams.GPURealTimeDisplayCoefficient*mean(abs(GPUParams.IQData),3));
    end

% end
return
end %-EF#5

function saveRFfastData_Information() 
%-EF#11

    % keyboard;
    tic;

    %%% Function ton only save the information of the acquisition
    
    % UserSet information
    UserSet=evalin('base','UserSet');
        % Voltage use during the acquisition
        hv1Sldr = findobj('Tag','hv1Sldr');
        UserSet.RealTime.Voltage=get(hv1Sldr,'Value');
        hv2Sldr = findobj('Tag','hv2Sldr');
        UserSet.THI.Voltage=get(hv2Sldr,'Value');
        clearvars hv1Sldr hv2Sldr;

        % Acquisition information
    Trans = evalin('base','Trans');
    P = evalin('base','P');
    Receive = evalin('base','Receive');
        Receive = Receive([Receive.bufnum] == 1);
        Receive = Receive(1);
    TX = evalin('base','TX');
    TW = evalin('base','TW');
    TPC = evalin('base','TPC');
    GPUParams = evalin('base','GPUParams');

    %%% Extract TimeFileInformation
    try 
        TimeFileInformation=evalin('base','TimeFileInformation');
    catch
        %%% we do not have information from motor
        disp('TimeFileInformation information set to 0');
        TimeFileInformation.info = 'NO INFORMATION';
    end

    %%% Extract TimeFileInformation
    try 
        Tracking=evalin('base','Tracking');
    catch
        %%% we do not have information from motor
        disp('Tracking information set to 0');
        Tracking.info = 'NO Tracking';
    end

    % Folder where to save generated from the acquisition button
    savedir = evalin('base','folder_acqui');
        RcvLastFrame=UserSet.RealTime.NumFrame;
    fname = evalin('base','filename_acqui');

    tic
    save([savedir,fname 'Information'],'UserSet','RcvLastFrame','Trans','P','Receive','TX','TW','TPC','GPUParams','TimeFileInformation','Tracking','-v7.3');
    t = toc;
    disp(['INFO saving in ',num2str(t),'s']);

    return
end %-EF#11
