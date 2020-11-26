function hypoglycemicEvents = findHypoglycemicEvents(data)
%findHypoglycemicEvents function that finds the hypoglycemic events in a 
%given glucose trace. The definition of hypoglycemic event can be found 
%in Danne et al.
%
%Input:
%   - data: a timetable with column `Time` and `glucose` containing the 
%   glucose data to analyze (in mg/dl).
%Output:
%   - hypoglycemicEvents: a structure containing the information on the
%   hypoglycemic events found by the function with fields:
%       - time: a vector containing the starting timestamps of each found 
%       hypoglycemic event;
%       - duration: a vector containing the duration (in min) of each found 
%       hypoglycemic event;
%
%Preconditions:
%   - data must be a timetable having an homogeneous time grid.
%
% ---------------------------------------------------------------------
%
% REFERENCE:
%   - Danne et al., "International consensus on use of continuous glucose
%   monitoring", Diabetes Care, 2017, vol. 40, pp. 1631-1640. DOI: 
%   10.2337/dc17-1600.
%
% ---------------------------------------------------------------------
%
% Copyright (C) 2020 Giacomo Cappon
%
% This file is part of AGATA.
%
% ---------------------------------------------------------------------
    
    %Check preconditions 
    if(~istimetable(data))
        error('findHypoglycemicEvents: data must be a timetable.');
    end
    if(var(seconds(diff(data.Time))) > 0)
        error('findHypoglycemicEvents: data must have a homogeneous time grid.')
    end
    
    
    k = 1; %hypoglycemicEvent vector current index
    nSamples = round(15/minutes(data.Time(2) - data.Time(1))); %number of consecutive samples required to define a valid event
    TH = 70; %set the hypoglycemic threshold 
    
    count = 0; %number of current found consecutive samples
    tempStartTime = []; %preallocate the hypo event starting time;
    hypoglycemicEvents.time = datetime(); %preallocate the result structure;
    hypoglycemicEvents.duration = []; %preallocate the result structure;
    
    flag = 0; %state flag -> 0: no event, 1: found a valid hypo event and 
              %currently in hypo, -1: found a valid hypo event and
              %currently not in hypo.
    
    for t = 1:height(data)
        
        if( data.glucose(t) < TH )
            
            %If it is a new event, reset count and set the hypothetical
            %starting time to the current timestamp
            if(count <= 0 )
                count = 0;
                tempStartTime = data.Time(t);
            end
            
            count = min([nSamples count + 1]); %limit count to nSamples
            
            %if count touches the "nSamples sample goal" we found a new event
            if(count == nSamples || flag == -1)
                flag = 1;
            end
            
        else
            
            
            if(flag == 0 && count > 0) 
                count = 0;
            end
            
            if(flag == 1)
                count = nSamples;
                flag = -1;
            end
            %count down
            
            count = count - 1;
            
        end
                
        if(count == 0 && flag == -1)
            hypoglycemicEvents.time(k) = tempStartTime;
            hypoglycemicEvents.duration(k) = minutes(data.Time(t-(nSamples-1)) - tempStartTime);
            k = k + 1;
            flag = 0;
        end

    end
    
    if(k == 1)
        hypoglycemicEvents.time = [];
    end
    
end
















