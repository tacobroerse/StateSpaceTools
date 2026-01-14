function[time,t,y,ntimes,syn]=synthetictimeseries(syn,synvar,settings)
%----------------------------------------------------------------------------   
% create synthetic time series
% settings contains settings
% synvar contains process variances for creating time series
%
%
% d.b.t.broerse@tudelft.nl
%
%----------------------------------------------------------------------------

    % some time series parameters
    ntimes=400;
    Step=100;
    time=[];
    % start
    
    if settings.continuoustime
        time(1)=0;
        for i=2:ntimes
            time(i)=time(i-1)+5/(ntimes)*(randn(1))^2;
        end
        time=time';
    else
        time=[1:ntimes]'/52;
    end
    
    Period=(time(end)-time(1))/(ntimes-1);
    t=(time-time(1))/Period+1;
    if ~isempty(settings.periodscycle)
        syn.lambda=2*pi/settings.periodscycle*Period;
    else
        %syn.lambda=zeros(size(t));
        syn.lambda=0;
    end
    if settings.continuoustime
        dt=diff(t);
        for i=1:ntimes
           y(i)=syn.mu(i)+syn.cycle(i,1)+sqrt(synvar.irr)*randn(1);
           if i<ntimes
               % local level
               syn.mu(i+1)=syn.mu(i)+syn.nu(i)*dt(i)+0.5*syn.xi(i)*dt(i)^2+sqrt(synvar.level)*dt(i)*randn(1);
               % slope
               syn.nu(i+1)=syn.nu(i)+syn.xi(i)*dt(i)+sqrt(synvar.slope)*randn(1)*dt(i);
               % acceleration
               syn.xi(i+1)=syn.xi(i)+sqrt(synvar.acc)*randn(1)*dt(i);
               if ~isempty(settings.periodscycle)
               syn.cycle(i+1,1)=[cos(syn.lambda*dt(i))  sin(syn.lambda*dt(i))]*syn.cycle(i,:)'+sqrt(synvar.cycle(1))*randn(1)*dt(i);
               syn.cycle(i+1,2)=[-sin(syn.lambda*dt(i)) cos(syn.lambda*dt(i))]*syn.cycle(i,:)'+sqrt(synvar.cycle(2))*randn(1)*dt(i); 
               else
syn.cycle(i+1,1)=0;
syn.cycle(i+1,2)=0;
               end
           end
        end
       % y(InterventionTime:end)=y(InterventionTime:end)+Step;
    
    else
        for i=1:ntimes
           y(i)=syn.mu(i)+syn.cycle(i,1)+sqrt(synvar.irr)*randn(1);         
           % local level
           syn.mu(i+1)=syn.mu(i)+syn.nu(i)+0.5*syn.xi(i)+sqrt(synvar.level)*randn(1);
           % slope
           syn.nu(i+1)=syn.nu(i)+syn.xi(i)+sqrt(synvar.slope)*randn(1);
           % acceleration
           syn.xi(i+1)=syn.xi(i)+sqrt(synvar.acc)*randn(1);
           syn.cycle(i+1,1)=[cos(syn.lambda)  sin(syn.lambda)]*syn.cycle(i,:)'+sqrt(synvar.cycle(1))*randn(1);
           syn.cycle(i+1,2)=[-sin(syn.lambda)  cos(syn.lambda)]*syn.cycle(i,:)'+sqrt(synvar.cycle(2))*randn(1);   
        end
    end
    figure; 
    hold on
    plot(time,y,'k-x')
    y=y';
    title('synthetic time series')
    
    if ~settings.continuoustime
        % Sampling frequency 
        Frequency=1/Period;
        % do Fourier analysis
        [f,Ypsd]=SpectralAnalysis(y,t,Frequency);
    else
        % lomb scargle
        [Ypsd,f]=lomb(y,time);
    end

    % Plot single-sided amplitude spectrum.
    figure
    semilogx(f,Ypsd,'k')
    %semilogx(f,Ypsd,'Color',Colors(ii,:)) 
    hold on
    title('Single-Sided Amplitude Spectrum of synthetic y(t)')
    xlabel('Frequency (1/yr)')
    ylabel('|Y(f)|')
end