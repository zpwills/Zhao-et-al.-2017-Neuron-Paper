function findCaPeaks(nSDs, numPeriods, deleteFocal)
%function [nPeaksRCaMP,nPeaksRhoA] = findCaPeak()
%This code assumes a particular directory structure, namely that the files
%are saved in directory Neuron1. You will need to modify the paths where
%the files are saved to get it working correctly.

%User defined now
%nSDs = 1; %number of standard deviations above mean to threshold by.
%numPeriods = 4; %number of time periods signal must exceed
%get user defined working directory' and cd to that directory
%get files matching pattern
folder_name = uigetdir('','Specify input directory');
folder_out_name = uigetdir('','Specify output directory');
files = dir([sprintf('%s/*.xlsx', folder_name)]);

for file = files'
    cd(folder_name);
    %Loop through sheet names?
    
    [status,sheets] = xlsfinfo(file.name);
    numOfSheets = numel(sheets); 
    [~,raw] = xlsread(file.name, 1); %gives col headers
    data = xlsread(file.name, 1); %gives raw data
    RCaMPcols = ~cellfun('isempty', strfind(raw, 'RCaMP')); %find RCaMP cols
    
    RhoAcols = ~cellfun('isempty', strfind(raw, 'RhoA')); %find RhoA cols
    RhoASensorNames = cellstr(raw(1,RhoAcols));
    RCaMPSensorNames = cellstr(raw(1,RCaMPcols));

    RCaMP = data(:,RCaMPcols);
    RhoA = data(:,RhoAcols);
    absTime = data(:,1);
    eventTimes = unique(xlsread(file.name, 2)*  86400); %% delete every second row
    

    %setting the Matlab figure
    f=figure('visible','off')

    %loop through all eventTimes, remove +/- 5 s around where user was focusing
    if deleteFocal;    
        for event = 1:length(eventTimes);
            filterFocus(:,event) = (data(:,1)<eventTimes(event)+5) & data(:,1)>(eventTimes(event)-5);
        end

        %instead of filtering out data,set it NaN 
        toFilter = sum(filterFocus,2);
        texpend = length(filterFocus);
        toFilter([texpend-5:length(filterFocus)],1) = 1;
        toFilter(toFilter>0 )=1;
        idx = find(toFilter==1);
        RCaMP_filtered(idx,:) = NaN;
        RhoA_filtered(idx,:) = NaN;
        data_filtered(idx,:) = NaN;
    end;


    %Demean Data (to do: high pass filter?)
    %deMeanedRhoA = bsxfun(@rdivide,RhoA,nanmean(RhoA));
    deMeanedRhoA = RhoA;
    deMeanedRCaMP = bsxfun(@rdivide,RCaMP,nanmean(RCaMP));
    sdRhoA = nanstd(deMeanedRhoA);
    sdRCaMP = nanstd(deMeanedRCaMP);
    outputData = struct;
    for sensor = 1:size(RCaMP,2);
        sensorData = deMeanedRCaMP(:,sensor); % grab the sensor
        n = 400; % average every n values

        % arbitrary data
        b = arrayfun(@(i) nanmean(sensorData(i:i+n-1)),1:n:length(sensorData)-n+1)'; % the averaged vector
        %findpeaks(sensorData, 'MinPeakProminence',.35)
        %pause
        %means = [];
        %for i = 1:length(b);
        %means = [means; repmat(b(i),n,1)];
        %end;
        %means = [means; repmat(means(end),length(sensorData)-length(means),1)];
        %means = medfilt1(sensorData,1);
        %means = repmat(nanmean(sensorData([1:endmeanCalc],1)), 1,length(means))';
        means = sgolayfilt(sensorData,1,121);
        %movavg = tsmovavg(sensorData,'s',100,1);
        tsig = sensorData >= (nSDs)*sdRCaMP(sensor) + means;
        %mean(sensorData([1:endmeanCalc],1))

        dsig = diff(tsig);
        startIndex = find(dsig > 0);
        endIndex = find(dsig < 0)-1;
        dsig(dsig==0) = [];

        %check if dsig is nonempty
        if length(dsig)>0;
            if dsig(1) == -1;%starts on fall.
                startIndex(2:length(startIndex)+1,:) =  startIndex(1:length(startIndex),:);
                startIndex(1) = 0;
            end
            if dsig(end) == 1;%ends on rise.
                endIndex(end+1) = length(sensorData);
            end
        end

        lenPeaks = (endIndex-startIndex)+1;
        supraPeaks = (lenPeaks >= numPeriods);

        peakStart =  startIndex(supraPeaks);%in absolute time
        peakEnd = endIndex(supraPeaks);
        for peak = 1: length(peakStart);

            startTime = peakStart(peak);
            endTime = peakEnd(peak);
            peakTimes = startTime:1:endTime;
            peakMags = sensorData([startTime:endTime],1);

            %Exclude if user requested of exclusion of peaks around focal events 
            if (deleteFocal & any(ismember(peakTimes, idx)));
                peakData{peak,1} = NaN;
                peakData{peak,2} = NaN;
            else    
                peakData{peak,1} = absTime(peakTimes);
                peakData{peak,2} = peakMags;

            end
        end;
        peakData_copy = peakData;
        peakData(cellfun(@(x) any(isnan(x)),peakData)) = {''}; % matlab is strange so need both these two lines. understandbly awkward
        peakData(all(cellfun(@isempty,peakData),2),:) = [];
        nPeaksRCaMP(:,sensor) = length(peakData);
        excludedPeaks(:,sensor) = length(peakData_copy) - length(peakData);
        roiName = sprintf('sensor%s', int2str(sensor));
        outputData.(roiName) = peakData;

        %Make figure for each sensor

        subplot(size(RCaMP,2),1, sensor); 
        plot(absTime,sensorData);
        %hline = refline([0 (nSDs)*sdRCaMP(sensor) + mean(sensorData(sensor))]);
        %set(hline,'LineStyle','--');
        hold on
        peakAbsTimes = cellfun(@(v) v(1), peakData(:,1));
        peakAbsMags = cellfun(@(v) v(1), peakData(:,2))
        %plot(peakAbsTimes, repmat((nSDs*sdRCaMP(sensor) + mean(sensorData(sensor))),length(peakAbsTimes),1), 'co', 'MarkerSize', 1.3)
        plot(peakAbsTimes, peakAbsMags, 'co', 'MarkerSize', 1.3)
        plot(absTime,(nSDs)*sdRCaMP(sensor) + means, '--');
        ylabel(sensor);
        set(gca,'FontSize',3) 
        axis([0,max(data(:,1)) + 10,0,4]);
        %title(sprintf('RCaMP Sensor %d', sensor));
        %ylabel('Delta F / F');
        %xlabel('Absolute Time (s)');
        clear peakData
    end
    fig = gcf
    paperUnits = get(fig, 'PaperUnits');
    set(fig,'PaperUnits','inches');
    paperSize = get(fig,'PaperSize');
    paperPosition = [.5 .5 paperSize - .5];
    set(fig,'PaperPosition', paperPosition);
    set(fig,'PaperUnits',paperUnits);
    % set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
    [~,filename,~] = fileparts(file.name);
    cd(folder_out_name);
    save(filename, 'outputData');

    figure_title = sprintf('%s.eps', filename);
    saveas(gcf,[figure_title], 'epsc' )
    close all    

    %Write data to directory:
    %create the directory dynamically. check if ~ exists, otherwise create
    %loop through all files in datain dir.

    nPeaksRCaMP = nPeaksRCaMP';
    excludedPeaks = excludedPeaks';

    TRCaMP = table(nPeaksRCaMP,excludedPeaks ,'RowNames', RCaMPSensorNames);
    writetable(TRCaMP, sprintf('%s_summary.csv', filename),'WriteRowNames',true )
    clearvars -except deleteFocal nSDs numPeriods endmeanCalc folder_name folder_out_name
    close all
end;


regionprops(BW,'centroid');
plot(ans.Centroid)
s = regionprops(BW,'centroid');
centroids = cat(1, s.Centroid);
imshow(BW)
hold on
plot(centroids(:,1),centroids(:,2), 'b*')