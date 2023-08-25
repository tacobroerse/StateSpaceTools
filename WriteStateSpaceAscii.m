function WriteStateSpaceAscii(settings,tseries,model,fitsmooth,estsmooth,varsmooth,distvar,nsample)
% save state space results as ascii

minarg=7;
maxarg=8;
narginchk(minarg,maxarg);

if nargin == minarg
    nsample = [];
end



if ~isfield(tseries,'tseriesname')
    tseries.tseriesname=tseries.name;
end


% make file for each variable in a multivariate model
for p=1:tseries.ntseries
    % formats
    clear headerstr formatstr1 formatstr2 writearray
    
    formatstr1='%15s';
    formatstr2='%15.8e';
    
    
    % header
    headerstr{1} = 'time [yr]';
    % data: time
    writearray = [tseries.time];
    
    % save files
    if isfield(tseries,'name')
        if tseries.ntseries==1
            savefiletimevariable=strcat(settings.maindir,'/',settings.savedir,'/statespaceresults_timevariable_',tseries.name,'.txt');
        else
            savefiletimevariable=strcat(settings.maindir,'/',settings.savedir,'/statespaceresults_timevariable_',tseries.name{p},'.txt');
        end
    else
        savefiletimevariable=strcat(settings.maindir,'/',settings.savedir,'/statespaceresults_timevariable_series_',num2str(p),'_sample_',num2str(nsample),'.txt');
    end
    % update header and array
    if tseries.ntseries==1
        
        headerstr=[headerstr,{tseries.tseriesname,'fit','residual','trend'}];
    else
        headerstr=[headerstr,{tseries.tseriesname{p},'fit','residual','trend'}];
    end
    writearray = [writearray,tseries.Y(p,:)',fitsmooth.fit{p}',estsmooth.epsilon(p,:)',fitsmooth.trend{p}'];
    formatstr1=strcat(formatstr1,' %15s %15s %15s  %15s');
    formatstr2=strcat(formatstr2,' %15.8e %15.8e %15.8e %15.8e');
    
    % add slope
    if settings.slope
        ii=length(headerstr)+1;
        headerstr{ii}='slope';
        writearray=[writearray fitsmooth.slope{p}'];
        
        formatstr1=strcat(formatstr1,' %15s');
        formatstr2=strcat(formatstr2,' %15.8e');
        
    end
    
    % add cycle standard
    if settings.cycle
        for i=1:settings.numbercycles
            ii=length(headerstr)+1;
            formatstr1=strcat(formatstr1,' %15s');
            formatstr2=strcat(formatstr2,' %15.8e');
            headerstr{ii}=strcat('cycle: ',num2str(settings.periodscycle(i)));
            writearray=[writearray fitsmooth.cycle{p}(i,:)'];
        end
    end
    
    if settings.AR
        ii=length(headerstr)+1;
        formatstr1=strcat(formatstr1,' %15s');
        formatstr2=strcat(formatstr2,' %15.8e');
        headerstr{ii}=strcat('AR');
        writearray=[writearray fitsmooth.AR{p}'];
    end
    
    % add trend standard deviation error
    ii=length(headerstr)+1;
    headerstr{ii}='trend stdev';
    stdevlevel=sqrt(varsmooth.level{p});
    writearray=[writearray stdevlevel'];
    formatstr1=strcat(formatstr1,' %15s');
    formatstr2=strcat(formatstr2,' %15.8e');
    
    if settings.slope
        ii=length(headerstr)+1;
        headerstr{ii}='slope stdev';
        stdevslope=sqrt(varsmooth.slope{p});
        writearray=[writearray stdevslope'];
        formatstr1=strcat(formatstr1,' %15s');
        formatstr2=strcat(formatstr2,' %15.8e');
    end
    
    % add cycle standard deviation error
    if settings.cycle
        for i=1:settings.numbercycles
            ii=length(headerstr)+1;
            
            formatstr1=strcat(formatstr1,' %15s');
            formatstr2=strcat(formatstr2,' %15.8e');
            headerstr{ii}=strcat('cycle: ',num2str(settings.periodscycle(i)),' stdev');
            stdevcycle=sqrt(squeeze(varsmooth.cycle{p}(1,:))');
            writearray=[writearray stdevcycle];
        end
    end

    % add AR standard deviation error
     if settings.AR
        ii=length(headerstr)+1;
        formatstr1=strcat(formatstr1,' %15s');
        formatstr2=strcat(formatstr2,' %15.8e');
        headerstr{ii}=strcat('AR stdev');
        stdevAR=sqrt(squeeze(varsmooth.AR{p})');
        writearray=[writearray stdevAR];
    end
    
    % end lines
    formatstr1=strcat(formatstr1,'\n');
    formatstr2=strcat(formatstr2,'\n');
    
    % print to file
    
    fid = fopen(savefiletimevariable,'w');
    fprintf(fid,formatstr1,headerstr{1:length(headerstr)});
    fprintf(fid,formatstr2,writearray');
    fclose(fid);
end



% write time invariant data
% save files

%if isfield(tseries,'name')
%    savefiletimeinvariable=strcat(settings.savedir,'/statespaceresults_disturbance_covariance_',tseries.name,'.txt');
%else
    savefiletimeinvariable=strcat(settings.maindir,'/',settings.savedir,'/statespaceresults_disturbance_covariances','.txt');
%end

% disturbance covariances

% distvar.level
clear headerstr
formatstr1='%15s';
formatstr2='%15.6e';
[ndim,~]=size(distvar.varcoveta);
writearray=[distvar.varcoveta];



% adjust format to size of covariance matrix
formatstr1=repmat([formatstr1,' '],1,ndim);
formatstr2=repmat([formatstr2,' '],1,ndim);


% make headers
for p=1:tseries.ntseries
    headerstr{model.indexlevel{p}}='level';
    if tseries.ntseries==1
        headertop{model.indexlevel{p}}=tseries.tseriesname;
    else
        headertop{model.indexlevel{p}}=tseries.tseriesname{p};
    end
    
    if settings.slope
        headerstr{model.indexslope{p}}=['slope '];
        if tseries.ntseries==1
            headertop{model.indexslope{p}}=tseries.tseriesname;
        else
            headertop{model.indexslope{p}}=tseries.tseriesname{p};
        end
    end
    if settings.acc
        headerstr{model.indexacc{p}}=['acc '];
        if tseries.ntseries==1
            headertop{model.indexacc{p}}=tseries.tseriesname;
        else
            headertop{model.indexacc{p}}=tseries.tseriesname{p};
        end
    end
    if settings.AR
        headerstr{model.indexAR{p}}=['AR '];
        if tseries.ntseries==1
            headertop{model.indexAR{p}}=tseries.tseriesname;
        else
            headertop{model.indexAR{p}}=tseries.tseriesname{p};
        end
    end
    if settings.cycle
        for i=1:settings.numbercycles
            
            headerstr{model.indexcycles{p}(i*2-1)}=['cycle per. ',num2str(settings.periodscycle(i))];
            headerstr{model.indexcycles{p}(i*2)}=['cycle* per. ',num2str(settings.periodscycle(i))];
            if tseries.ntseries==1
                headertop{model.indexcycles{p}(i*2-1)}=tseries.tseriesname;
                headertop{model.indexcycles{p}(i*2)}=tseries.tseriesname;
            else
                headertop{model.indexcycles{p}(i*2-1)}=tseries.tseriesname{p};
                headertop{model.indexcycles{p}(i*2)}=tseries.tseriesname{p};
            end
        end
    end
    
end

% end lines
formatstr1=strcat(formatstr1,'\n');
formatstr2=strcat(formatstr2,'\n');

% print to file

fid = fopen(savefiletimeinvariable,'w');
fprintf(fid,formatstr1,headertop{1:length(headertop)});
fprintf(fid,formatstr1,headerstr{1:length(headerstr)});
fprintf(fid,formatstr2,writearray');
fclose(fid);



%% irregular variance
if isfield(tseries,'name')
    savefiletimeinvariable2=strcat(settings.maindir,'/',settings.savedir,'/statespaceresults_irregular_covariance_',tseries.generalname,'.txt');
else
    savefiletimeinvariable2=strcat(settings.maindir,'/',settings.savedir,'/statespaceresults_irregular_covariance_',num2str(nsample),'.txt');
end

% disturbance covariances

% distvar.level
clear headerstr
formatstr1='%15s';
formatstr2='%15.6e';

writearray=[distvar.varcovirr];



% adjust format to size of covariance matrix
formatstr1=repmat([formatstr1,' '],1,tseries.ntseries);
formatstr2=repmat([formatstr2,' '],1,tseries.ntseries);


% make headers
for p=1:tseries.ntseries
    if tseries.ntseries == 1
        headertop2=tseries.tseriesname;
    else
        headertop2=tseries.generalname;
    end
    
end

% end lines
formatstr1=strcat(formatstr1,'\n');
formatstr2=strcat(formatstr2,'\n');

% print to file

fid = fopen(savefiletimeinvariable2,'w');
fprintf(fid,formatstr1,headertop2(1:length(headertop2)));
fprintf(fid,formatstr2,writearray');
fclose(fid);


% % save time invariant: distvar.slope, distvar.irregular,
% % distvarcycle, trend, trend_sigma, LogL, convergence error
% % formats
% formatstr1='%30s %30s %30s %30s %30s %30s';
% formatstr2='%30.8e %30s %30.8e %30.8e %30.8e %30.8e';
%
% % header
% headerstr = {'# log likelihood', 'convergence error', 'mean trend', 'mean trend stdev', 'slope disturbance variance','irregular variance'};
% slopestdev=sqrt(meanslope.slopevar);
% writecell = {estfilter.LogL,ConvError.type, meanslope.slope , slopestdev, distvar.slope, distvar.irr};
% % add cycle variance
% if settings.cycle
%     for i=1:settings.numbercycles
%         % 1st component
%         ii=length(headerstr)+1;
%         headerstr{ii}=strcat('cycle: ',num2str(settings.periodscycle(i)),' disturbance variance ');
%
%         writecell{ii}=distvar.cycle(1,i);
%         % 2nd component: *
%         ii=length(headerstr)+1;
%         headerstr{ii}=strcat('cycle: ',num2str(settings.periodscycle(i)),' disturbance variance *');
%         writecell{ii}=distvar.cycle(2,i);
%         formatstr1=strcat(formatstr1,' %30s %30s');
%         formatstr2=strcat(formatstr2,' %30.8e %30.8e');
%     end
% end
%
%
%
% % end lines
% formatstr1=strcat(formatstr1,'\n');
% formatstr2=strcat(formatstr2,'\n');
%
% % print to file
%
% fid = fopen(savefiletimeinvariable,'w');
% fprintf(fid,formatstr1,headerstr{1:length(headerstr)});
% fprintf(fid,formatstr2,writecell{1:length(writecell)});
% fclose(fid)

end

