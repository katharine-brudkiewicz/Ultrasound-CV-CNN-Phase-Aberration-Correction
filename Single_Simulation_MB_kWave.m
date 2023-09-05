%% k-Wave Simulation of Transcranial ULM with Micrbubbles Simulated with Single Simulation Method

% requires k-Wave and USTB toolboxes to be installed


clearvars
close all;
for k=1:215
    tic
    clearvars -except k
    close all;
    %% Basic definitions
    %
    % Constants
    
    f0 = 2e6;       % pulse center frequency [Hz]
    cycles=2;       % number of cycles in pulse
    c0 = 1540;      % medium speed of sound [m/s]
    rho0 = 1020;    % medium density [kg/m3]
    F_number = 0.7; % F number for CPWC sequence (i.e. maximum angle)
    N=1;            % number of plane waves in CPWC sequence
    
    %
    % Linear Array Probe defined as a USTB structure.
    
    prb=linear_array();
    prb.N=128;                   % number of elements 
    prb.pitch=300e-6;           % probe pitch in azimuth [m]
    prb.element_width=270e-6;   % element width [m]
    prb.element_height=4.5e-3; % element height [m]
    
    fig_handle = prb.plot([],'Linear array');
    
    %% Computational grid
    %
    f_max = 1.2*f0;
    lambda_min = c0/f_max;
    
    % mesh resolution, choose one
    mesh_resolution='element4'; 
    switch mesh_resolution
        case 'element2' % around 8 min per wave
            dx=prb.pitch/2;                                         % 2 elements per pitch 
        case 'element4' % around 52 min per wave
            dx=prb.pitch/3;                                         % 2 elements per pitch 
        otherwise
            error('Not a valid option');
    end
    
    % mesh size
    PML_size = 20;                                          % size of the PML in grid points
    Nx=round(40e-3/dx); Nx=Nx+mod(Nx,2);
    Nz=round(40e-3/dx); Nz=Nz+mod(Nz,2);
    grid_width=Nx*dx;
    grid_depth=Nz*dx;
    domain=linear_scan('x_axis', linspace(-grid_width/2,grid_width/2,Nx).', 'z_axis', linspace(0,grid_depth,Nz).');
    
    kgrid = kWaveGrid(domain.N_z_axis, domain.z_step, domain.N_x_axis, domain.x_step);
    
    %% Propagation medium
    %
    scattering_map = randn([Nx, Nz]);
    %% Bubble Speed of Sound/Density
    cbubble = 500
    rhobubble = 10
    scattering_c0 = cbubble* scattering_map;
    scattering_c0(scattering_c0 > 650) = 650;
    scattering_c0(scattering_c0 < 350) = 350;
    scattering_rho0 = scattering_c0 / 50;
    
    % Skull speed of sound and density
    cskullcortical  = 2800;     % speed of sound [m/s]
    rhoskullcortical= 1850;
    cskulltrabecular= 2300; % speed of sound [m/s]
    rhoskulltrabecular= 1700; 
    alphabrain = 1.2;
    alphatrab = 32;
    alphacort = 16;
    alpha0 = 0;
    
   
    scattering_cb = cskullcortical*scattering_map;
    scattering_cb(scattering_cb > 3000) = 3000;
    scattering_cb(scattering_cb < 2600) = 2600;
    scattering_rhob=scattering_cb/1.5;
    
% Import Skull segment    
    
    patientIDs = 'patient';
    load([patientIDs num2str(k,'%d') '.mat'],'skull');
    ttt = randi([20 290])
    segment = skull(:,:,ttt)
    segment1 = rot90(segment)
    mediummm = flip(segment1,1);
    imagesc(mediummm)
    axis equal
    colorbar
    title('Acoustic Pressure field source')
    xlabel('Axial Position [mm]');
    ylabel('Lateral Position [mm]');

    % define properties
    sound_speed_map = c0 * ones(Nx, Nz); %.* background_map;
    density_map = rho0 * ones(Nx, Nz); %.* background_map;

    background_map_mean = 1;
    background_map_std = 0.004;
    background_map = background_map_mean + background_map_std * randn([Nx, Nz]);
    

    sound_speed_map = c0 * ones(Nx, Nz) .* background_map;
    density_map = rho0 * ones(Nx, Nz) .* background_map;
    
%     %define bubbles for a highly scattering region
    radius = 5e-5;
    y_pos = 8e-3;
    x_pos = 30e-3;
    scattering_region1 = makeDisc(Nx, Nz, round(x_pos / dx), round(y_pos / dx), round(radius/dx));
    
    y_pos = 8e-3;
    x_pos = 35e-3;
    scattering_region2 = makeDisc(Nx, Nz, round(x_pos / dx), round(y_pos / dx), round(radius/dx));
    
    
    y_pos = 16e-3;
    x_pos = 30e-3;
    scattering_region3 = makeDisc(Nx, Nz, round(x_pos / dx), round(y_pos / dx), round(radius/dx));
    
    y_pos = 16e-3;
    x_pos = 35e-3;
    scattering_region4 = makeDisc(Nx, Nz, round(x_pos / dx), round(y_pos / dx), round(radius/dx));
    
    y_pos = 24e-3;
    x_pos = 30e-3;
    scattering_region5 = makeDisc(Nx, Nz, round(x_pos / dx), round(y_pos / dx), round(radius/dx));
    % %     
    y_pos = 24e-3;
    x_pos = 35e-3;
    scattering_region6 = makeDisc(Nx, Nz, round(x_pos / dx), round(y_pos / dx), round(radius/dx));
    
    y_pos = 32e-3;
    x_pos = 30e-3;
    scattering_region7 = makeDisc(Nx, Nz, round(x_pos / dx), round(y_pos / dx), round(radius/dx));
    
    y_pos = 32e-3;
    x_pos = 35e-3;
    scattering_region8 = makeDisc(Nx, Nz, round(x_pos / dx), round(y_pos / dx), round(radius/dx));
    %
    y_pos = 20e-3;
    x_pos = 25e-3;
    scattering_region9 = makeDisc(Nx, Nz, round(x_pos / dx), round(y_pos / dx), round(radius/dx));
    %     
    %     
    %assign region
    sound_speed_map(scattering_region1 == 1) = scattering_c0(scattering_region1 == 1);
    density_map(scattering_region1 == 1) = scattering_rho0(scattering_region1 == 1);
% %     
% %  assign bubbles to medium

    sound_speed_map(scattering_region2 == 1) = scattering_c0(scattering_region2 == 1);
    density_map(scattering_region2 == 1) = scattering_rho0(scattering_region2 == 1);
    sound_speed_map(scattering_region3 == 1) = scattering_c0(scattering_region3 == 1);
    density_map(scattering_region3 == 1) = scattering_rho0(scattering_region3 == 1);
    sound_speed_map(scattering_region4 == 1) = scattering_c0(scattering_region4 == 1);
    density_map(scattering_region4 == 1) = scattering_rho0(scattering_region4 == 1);
    sound_speed_map(scattering_region5 == 1) = scattering_c0(scattering_region5 == 1);
    density_map(scattering_region5 == 1) = scattering_rho0(scattering_region5 == 1);
    sound_speed_map(scattering_region6 == 1) = scattering_c0(scattering_region6 == 1);
    density_map(scattering_region6 == 1) = scattering_rho0(scattering_region6 == 1);
    sound_speed_map(scattering_region7 == 1) = scattering_c0(scattering_region7 == 1);
    density_map(scattering_region7 == 1) = scattering_rho0(scattering_region7 == 1);
     sound_speed_map(scattering_region8 == 1) = scattering_c0(scattering_region8 == 1);
     density_map(scattering_region8 == 1) = scattering_rho0(scattering_region8 == 1);
    %   
    sound_speed_map(scattering_region9 == 1) = scattering_c0(scattering_region8 == 1);
    density_map(scattering_region9 == 1) = scattering_rho0(scattering_region9 == 1);
    %   
    sound_speed_map(mediummm == 1) = scattering_cb(mediummm == 1);
    density_map(mediummm == 1) = scattering_rhob(mediummm == 1);
    %     
    % assign to the medium inputs
    medium.sound_speed = sound_speed_map;
    medium.density = density_map;
    % 

    %medium.alpha_coeff = 2;
    medium.alpha_power = 2; % default
    medium.alpha_coeff = 0.6 * ones(Nx, Nz); % general tissue

    figure;
    subplot(1,2,1);
    imagesc(domain.x_axis*1e3,domain.z_axis*1e3,medium.sound_speed); colormap gray; colorbar; axis equal tight;
    xlabel('x [mm]');
    ylabel('z [mm]');
    title('c_0 [m/s]');
    subplot(1,2,2);
    imagesc(domain.x_axis*1e3,domain.z_axis*1e3,medium.density); colormap gray; colorbar; axis equal tight;
    xlabel('x [mm]');
    ylabel('z [mm]');
    title('\rho [kg/m^3]');
    
    %% Time vector
    %
    % We define the time vector depending on the CFL number, the size of the
    % domain and the mean speed of sound.
    
    cfl=0.3;
    t_end=2*sqrt(grid_depth.^2+grid_depth.^2)/mean(medium.sound_speed(:));
    kgrid.makeTime(medium.sound_speed,cfl,t_end);
    
    %% Sequence
    %
    % We define a sequence of plane-waves
    alpha_max=1/2/F_number;                         % maximum angle span [rad]
    if N>1
        angles=linspace(-alpha_max,alpha_max,N);    % angle vector [rad]
    else
        angles = 0;
    end
    seq=wave();
    for n=1:N
        seq(n)=wave();
        seq(n).apodization = apodization('f_number',1,'window',window.rectangular,'focus',scan('xyz',[0 0 10e-3]));
        seq(n).source.azimuth=angles(n);
        seq(n).source.distance=-Inf;
        seq(n).probe=prb;
        seq(n).sound_speed=1540;    % reference speed of sound [m/s]
        seq(n).delay = min(seq(n).delay_values);
        seq(n).source.plot(fig_handle);
    end
    
    %% Source & sensor mask
    %
    % Based on the uff.probe we find the pixels in the domain that must work as
    % source and sensors.
    
    % find the grid-points that match the element
    source_pixels={};
    element_sensor_index = {};
    n=1;
    for m=1:prb.N_elements
        plot((prb.x(m)+[-prb.width(m)/2 prb.width(m)/2])*1e3,[0 0],'k+-'); hold on; grid on;
        source_pixels{m}=find(abs(domain.x-prb.x(m))<prb.width(m)/2 & abs(domain.y-prb.y(m))<prb.height(m) & abs(domain.z-prb.z(m))<=domain.z_step/2);
        element_sensor_index{m} = n:n+numel(source_pixels{m})-1;
        n=n+numel(source_pixels{m});
    end
    
    % sensor mask
    sensor.mask = zeros(domain.N_z_axis, domain.N_x_axis);
    for m=1:prb.N_elements
        sensor.mask(source_pixels{m}) = sensor.mask(source_pixels{m}) + 1;
    end
    
    % source mask
    source.u_mask=sensor.mask;
    
    figure;
    h=pcolor(domain.x_axis,domain.z_axis,source.u_mask); axis equal tight;
    title('Source/Sensor mask')
    set(h,'edgecolor','none');
    set(gca,'YDir','reverse');
    xlabel('x [mm]');
    ylabel('z [mm]');
    
    %% Calculation
    %
    % We are ready to launch the k-Wave calculation
    
    disp('Launching kWave. This can take a while.');
    for n=1:N
        delay=seq(n).delay_values-seq(n).delay;
        denay=round(delay/kgrid.dt);
        seq(n).delay = seq(n).delay - cycles/f0/2;
        
        % offsets
        tone_burst_offset = [];
        for m=1:prb.N_elements
            tone_burst_offset = [tone_burst_offset repmat(denay(m),1,numel(source_pixels{m}))];
        end
        source.ux = toneBurst(1/kgrid.dt, f0, cycles, 'SignalOffset', tone_burst_offset);   % create the tone burst signals
        source.uy = 0.*source.ux;
        source.u_mode ='dirichlet';
        
        % set the input arguements: force the PML to be outside the computational
        % grid; switch off p0 smoothing within kspaceFirstOrder2D
        input_args = {'PMLInside', false, 'PMLSize', PML_size, 'PlotPML', false, 'Smooth', false};
        
        % run the simulation
        sensor_data(:,:,n) = permute(kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:}),[2 1]);
    end
    

    sensor_data(isnan(sensor_data))=0;
    
    %% Gather element signals
    %
    % After calculaton we combine the signal recorded by the sensors according to the
    % corresponding element
    element_data=zeros(numel(kgrid.t_array),prb.N_elements,numel(seq));
    for m=1:prb.N_elements
        if  ~isempty(element_sensor_index{m})
            element_data(:,m,:)=bsxfun(@times,sqrt(1./kgrid.t_array).',trapz(kgrid.y(source_pixels{m}),sensor_data(:,element_sensor_index{m},:),2));
        end
    end
      
    %% Band-pass filter
    %
    % We remove some numerical noise by band-pass filtering
    filtered_element_data=band_pass(element_data,1/kgrid.dt,[0 1e6 8e6 10e6]);
    
    %% Channel_data
    %
    % We can now store the simulated data into a uff.channel_data class
    channel_data = channel_data();
    channel_data.probe = prb;
    channel_data.sequence = seq;
    channel_data.initial_time = 0;
    channel_data.sampling_frequency = 1/kgrid.dt;
    channel_data.data = filtered_element_data;
    
    % taking care of NaNs
    channel_data.data(isnan(channel_data.data))=0;
    
    %% Beamforming
    %
    % To beamform we define a new (coarser) uff.linear_scan. We also define the
    % processing pipeline and launch the beamformer
    
    domain=linear_scan('x_axis',linspace(domain.x_axis(1),domain.x_axis(end),512).',...
    'z_axis',linspace(domain.z_axis(1),domain.z_axis(end),512).');
    
    pipe = pipeline();
    pipe.channel_data = channel_data;
    pipe.scan = domain;
    
    pipe.receive_apodization.window = window.hanning;
    pipe.receive_apodization.f_number = F_number;
    
    b=coherent_compounding();
    b.dimension = dimension.both;
    
    das = pipe.go({das b});
    das.plot([],'DAS'); hold on;
    patientIDs = 'extra';
    
    dir1 = "D:\Ultrasound_Data\";
    m = k +645;
    filename=strcat(dir1,patientIDs,num2str(k,'%d'));
    save(filename)
    toc
end

