%% Marmottant simulations for Microbubble SonoVue Response

for k=1:215
    clearvars -except k
    close all;
    
    patientIDs = 'marmottantinput';
    
    dir1 = "D:\Input_Signal\";
    filename=strcat(dir1,patientIDs,num2str(k,'%d'));
    load(filename,"pressure","timeperiod")
    
    
    
    % define waveform parameters
    freq = 2;    % incident frequency in [MHz]
    ncyc = 3;    % number of cycles
    pnp = 10;   % desired peak negative pressure in [kPa]
    pol = 1;    % = 1 or -1, polarity of pulse
    
    % generate waveform
    [p0,p0t] = generate_waveform(freq,ncyc,pnp,pol);
    dT=p0t(2)-p0t(1); 
    f=(0:length(p0t)-1)/(length(p0t)-1)/dT;
    
    plot(f,dbzero(abs(fft(p0))))
    
    
    plot(p0t,p0);
    plot(timeperiod,pressure(1,:));
    
    %% marmottant simulation
    
    % define input parameters 
    distance = 50*1e-6;     % distance to compute scattered pressure in [m]
    bubbleType = 'literature';  % type of bubble to determine shell parameters
    R0 = 5;           % initial bubble radius in [um]
    freq = 2;            % incident frequency in [MHz]
    
    bubble_pressure = zeros(10,length(timeperiod));
    for g=1:10
        tic
        [time,rad,pscat,psurf] = marmottant(distance,bubbleType,R0,freq,pressure(g,:),timeperiod);
        toc
    
        %% interpolate for even time steps
        pscat2=interp1(time,pscat,timeperiod);
        dT=timeperiod(2)-timeperiod(1); 
        bubble_pressure(g,:) = pscat2;
    end
    
    patientIDs = 'marmottantresponse';
    
    dir1 = "D:\Marmottant_Response\";
    filename=strcat(dir1,patientIDs,num2str(k,'%d'));
    save(filename,"bubble_pressure")
end

