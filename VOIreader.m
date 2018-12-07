function [voi_points] = VOIreader(varargin)
%% PARSE INPUT
% in the future fix this to accept file and also ask for corresponding MRI file
if length(varargin) == 1
    filename = varargin{1};   
else
    filename = uigetfile('C:\','Select VOI file');
end

%% PARSE FILE

text = string(fileread(filename));
voi_text = strsplit(text,{'\n','\t'})';
voi_bool = 0;

sliceData = {}; %see comments below on how these are filled
%first column = slice #, second column = num points, third column = matrix of values

for i = 1:length(voi_text)
    rowText = voi_text(i);
    %check for VOI indicator
        if(length(regexp(rowText,'MIPAV VOI FILE','match'))>0)
            voi_bool = 1;
        end
    %find number of slices for VOI
        if(length(regexp(rowText,'number of slices for the VOI','match'))>0)
            totSlices_val = voi_text(i-1);
        end
    %grab data from slices
        if(length(regexp(rowText,'slice number','match'))>0)
            dummySlice = cell(1,3);
            %this is dependent on MIPAV format
            dummySlice{1} = str2num(voi_text(i-1)); %slice #
            dummySlice{2} = str2num(voi_text(i+3)); %num points
            dummySlice{3} = voi_text(i+5:i+5+str2num(voi_text(i+3))-1); %matrix of values
            sliceData = cat(1, sliceData, dummySlice);
        end    
end

%% CHECK OUTPUTS

%find 'MIPAV VOI FILE' otherwise throw error
    if(voi_bool == 0)
        disp('WARNING: FILE NOT MIPAV VOI FORMAT. EXITING.');
        return;
    end

%check to make sure the sliceData matches the same size as number of slices for the VOI
    if(str2num(totSlices_val) ~= size(sliceData,1))
        disp('WARNING: NUM SLICES DOES NOT MATCH EXTRACTED DATA');
        disp(['# slices: ' totSlices_val]);
        disp(['# slices with extracted data: ' int2str(size(sliceData,1))]);
    end

%% IF ALL CHECKS OUT -> CONVERT TO USEABLE FORMAT

% slice number in VOI starts at zero -- adjust to correct value
% convert string points to numeric and round
voi_points = cell(size(sliceData,1),2);

for i = 1:size(sliceData,1)
    voi_points{i,1} = sliceData{i,1}+1; %slice number
    sliceData_i = sliceData{i,3}; %grab string voxel data
    dummyPts = zeros(size(sliceData_i,1),2);       
    for j = 1:size(sliceData_i,1)
        points_j = regexp(sliceData_i(j,1),' ','split'); %split to x and y
        dummyPts(j,1) = round(str2num(points_j(1))); %x
        dummyPts(j,2) = round(str2num(points_j(2))); %y
    end
    voi_points{i,2} = dummyPts; %allocate points
end


