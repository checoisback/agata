function timeDelay = timeDelay(data,dataHat,ph)
%timeDelay function that computes the delay of a predicted glucose trace. 
%The time delay is computed as the time shift necessary to maximize the 
%correlation between the two traces.
%
%Inputs:
%   - data: a timetable with column `Time` and `glucose` containing the 
%   glucose data to analyze (in mg/dl);
%   - data: a timetable with column `Time` and `glucose` containing the inferred 
%   glucose data (in mg/dl) to compare with `data`;
%Output:
%   - timeDelay: the computed delay (min).
%
%Preconditions:
%   - data and dataHat must be a timetable having an homogeneous time grid;
%   - data and dataHat must contain a column named `Time` and another named `glucose`;
%   - data and dataHat must start from the same timestamp;
%   - data and dataHat must end with the same timestamp;
%   - data and dataHat must have the same length.
%
% ------------------------------------------------------------------------
% 
% Reference:
%   - Wikipedia on time delay using correlation. 
%   https://en.wikipedia.org/wiki/Cross-correlation#Time_delay_analysis (Accessed: 2020-12-10).
% 
% ------------------------------------------------------------------------
%
% Copyright (C) 2020 Giacomo Cappon
%
% This file is part of AGATA.
%
% ---------------------------------------------------------------------
    
    %Check preconditions 
    if(~istimetable(data))
        error('timeDelay: data must be a timetable.');
    end
    if(var(seconds(diff(data.Time))) > 0 || isnan(var(seconds(diff(data.Time)))))
        error('timeDelay: data must have a homogeneous time grid.')
    end
    if(~istimetable(data))
        error('timeDelay: dataHat must be a timetable.');
    end
    if(var(seconds(diff(data.Time))) > 0)
        error('timeDelay: dataHat must have a homogeneous time grid.')
    end
    if(data.Time(1) ~= dataHat.Time(1))
        error('timeDelay: data and dataHat must start from the same timestamp.')
    end
    if(data.Time(end) ~= dataHat.Time(end))
        error('timeDelay: data and dataHat must end with the same timestamp.')
    end
    if(height(data) ~= height(dataHat))
        error('timeDelay: data and dataHat must have the same length.')
    end
    if(~any(strcmp(fieldnames(data),'Time')))
        error('timeDelay: data must have a column named `Time`.')
    end
    if(~any(strcmp(fieldnames(data),'glucose')))
        error('timeDelay: data must have a column named `glucose`.')
    end
    if(~any(strcmp(fieldnames(dataHat),'Time')))
        error('timeDelay: dataHat must have a column named `Time`.')
    end
    if(~any(strcmp(fieldnames(dataHat),'glucose')))
        error('timeDelay: dataHat must have a column named `glucose`.')
    end
    
    %Get indices having no nans in both timetables
    idx = find(~isnan(dataHat.glucose) & ~isnan(data.glucose));
    
    % Compute metric
    timeDelay = finddelay( data.glucose(idx), dataHat.glucose(idx)) * minutes(data.Time(2)-data.Time(1));
    
    % remove missing data and find row indexes (TF)
    x = [data.glucose, dataHat.glucose];
    [~, TF] = rmmissing(x);

    % split data in regions when NaN (described by TF)
    idd = diff([true;TF;true]);
    idb = find(idd<0);
    ide = find(idd>0)-1;
    regions = arrayfun(@(b,e)x(b:e,:),idb,ide,'uni',0);

    count = 0;
    for r = 1:size(regions,1)
        
        if length(regions{r})>288
            
            count = count + 1;
            timeDelay(count,1) = getDelay(regions{r,1}(:,1),regions{r,1}(:,2),ph)*minutes(data.Time(2)-data.Time(1));
            
        end
    end
    
    timeDelay = mean(timeDelay);
    
end



function delay = getDelay(cgm_double,yhat_double,ph)

% delay

N=length(cgm_double);
% count=0;
for j = 0:ph
    
    for i = 1:N-ph-j
        indx = i+ph+j;
        dif(i) = (yhat_double(indx)-cgm_double(i+ph))^2;
    end
    count = nansum(dif);
    delay_tmp(j+1) = (1/(N-ph+1))*count;
    count = 0;
    
end

[~,d] = min(delay_tmp);
delay = d-1; % caused by we are using j+1 as index
end
