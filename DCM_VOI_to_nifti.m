% DCM_VOI_to_nifti.m
% requiresL: image processing toolbox
%
% Authors: Stephanie Harmon, Jonathan Sackett
% Institution: National Institutes of Health, National Cancer Institute, Molecular Imaging Branch
%              from the Lab of Dr. Ismail Baris Turkbey, MD
% Date: 10/8/2018
%
% given the path to a directory with proper file structure which contains CT data and VOI data
% will produce a nifti file from the region of interest defined by the VOI file.
% output image contains only the ROI plus a specified amount of padding. 
% output is placed in folder labeled nifti and output file is labeled with the name of the VOI file
% and the name of the parent directory.
%
% Arguments:
%   dcm_path = the path to the directory that contains the dicom files
%   voi_path = full path to voi file used for segmentation of your dicom files.
%   nifti_path = full pathname (including filename) to location that you want nifti file saved
%   padding = number of pixels outside the borders of the voi segmentation you would like to include


function DCM_VOI_to_nifti(dcm_path, voi_path, nifti_path, padding)

    %read in CT data
    dcm = readSeries(dcm_path); 
    % TODO: make this step take less time. 
    % this step takes the majority of the time to complete.
    %dcm = ConvertDicom(dcm_path); % other variant of a function to read in volumes
    %disp('Dicom converted')

    %read in VOI data
    voi = VOIreader(voi_path);
    %disp('VOI decoded')

    mask = create_mask(voi, size(dcm.data));
    %disp('mask created')

    shoebox = make_shoebox(mask); 
    %disp('shoebox completed')

    masked_image = apply_mask_shoebox(dcm, mask, shoebox, padding); %adds 5 pixels of padding to voi image.
    %disp('Masquerade!')

    niftiwrite(masked_image,nifti_path);
    %disp('Nifty!')

end

%###############################################################################

function mask = create_mask(voi,im_size)
    mask = zeros(im_size);
    for i = 1:size(voi,1)
        slice_pts = voi{i,2}; %fill slice points based on VOI file
        slice_z = voi{i,1};
        slice_mask = poly2mask(slice_pts(:,1),slice_pts(:,2),im_size(1),im_size(2)); %create mask from points
        mask(:,:,slice_z) = slice_mask; %fill mask at every slice
    end
end

%###############################################################################

function masked_image = apply_mask_shoebox(dcm, mask, sbx, pad)
% adds padding to mask and shoebox, then crops the image and the mask 
% according to the coordinates passed in 'sbx', then applies the mask to the 
% cropped image.

    shoebox{1} = [(sbx(1,1) - pad):(sbx(1,2) + pad)];
    shoebox{2} = [(sbx(2,1) - pad):(sbx(2,2) + pad)];
    shoebox{3} = [(sbx(3,1) - pad):(sbx(3,2) + pad)];

    cropped_image = dcm.data(shoebox{1},shoebox{2},shoebox{3});
    cropped_mask = mask(shoebox{1},shoebox{2},shoebox{3});
    
    se = strel('sphere',pad);
    cropped_mask = imdilate(cropped_mask, se);

    masked_image = cropped_mask .* cropped_image;
end

%###############################################################################

function  shoebox_faces = make_shoebox(mask)
% Given a 3D Binary mask will calculate the coordinates corresponding to the 
% the centers of the faces of a box with minimum size that would exactly contain
% the 3D mask.

    shoebox_faces = zeros(3,2);
    for i = 1:3 % loops once for each dimension
        a = size(mask);
        back_face = false;
        for j = 1:a(i) % loops once for each slice in each dimension
            switch i % selects which dimension to loop through
                case 1
                    slice_check = mask(j,:,:);
                case 2
                    slice_check = mask(:,j,:);
                case 3
                    slice_check = mask(:,:,j);
            end

            % records the 'first' slice containing a non-zero value
            if and(any(any(slice_check)),(back_face == false))
                shoebox_faces(i,1) = j;
                back_face = true;
            end
            % records the 'last' slice containing a non-zero value
            if and(all(all(not(slice_check))),(back_face == true))
                shoebox_faces(i,2) = (j-1);
                break;
            end
        end
    end
end

%###############################################################################

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
end

%###############################################################################
% TODO: Decide which function is faster. "readSeries" or "ConvertDicom". I suspect they are similar speeds. Unfortunately loops in matlab are slow but there is no way to vectorize a data import function. 

% FIXME: figure out why this exact same function runs so much faster in python, despite the fact that I've pre-allocated space for all variables in this version but not in the python version. 

function  out = readSeries(dirName)
% given the path to a directory containing DICOM files, reads in the DICOM files
% sorts the images in the DICOM file according the ImagePositionPatient[3]
% and returns a 3D matrix representing the voxel intensity at each coordinate.
        
        list = dir([dirName '\']);
        list = list(~strcmp({list.name},'.') & ~strcmp({list.name},'..'));
        list = {list(:).name};
        zpos = zeros(size(list));
        data = zeros([size(dicomread([dirName filesep list{1}])),numel(list)]);

        % loop takes ~0.08s per loop.
        % total time is 27.94s for a 512x512x332 file. 
        % 0.08s * 332 = 26.56s
        for i = 1:numel(list)
            file = list(i);
            filepath = [dirName filesep file{:}];
            data(:,:,i) = dicomread(filepath);
            info(i) = dicominfo(filepath);          
        end

        pos = horzcat(info(:).(dicomlookup('0020','0032'))).';
        zpos(:)=pos(:,3);
        [~,ordering] = sort(zpos,2);
        data_sorted = data(:, :, ordering);
        out = struct('data', data_sorted);

end

%###############################################################################

function [Geometry, header] = ConvertDicom(varargin)

    %% Parse Input %%
    
    if length(varargin) == 1
        directory = varargin{1};   
    else
        directory = uigetdir('C:\','Select directory with DICOM files');
    end
    
    %% READ IN DATA 
    
    %loop through each file in directory (ie, slice)
    %list = dir([directory '\*.dcm']);
    list = dir([directory '\']);
    list = list(~strcmp({list.name},'.') & ~strcmp({list.name},'..'));
    filenames = cellfun(@(x)fullfile(directory, x), {list.name}, 'uni', 0)';
    list = list(cellfun(@isdicom, filenames));
    
    info = cell(1,size(list,1));
    zpos = zeros(1,size(list,1)); %axial position, for sorting
    
    %read in all headers and z-positions
    for i = 1:size(list,1) %skip "." and ".."
        info{i} = dicominfo([directory '\' list(i).name]);
        zpos(i) = info{i}.ImagePositionPatient(3);
        ypos(i) = info{i}.ImagePositionPatient(2);
        xpos(i) = info{i}.ImagePositionPatient(1);
    end
    
    [B,ordering] = sort(zpos,2);
    
    %%
    
    %Number of slices, pixel size
    Nx = double(info{ordering(1)}.Width);
    Ny = double(info{ordering(1)}.Height);
    dx = info{ordering(1)}.PixelSpacing(1); %in mm
    dy = info{ordering(1)}.PixelSpacing(2); %in mm
    x_start = info{ordering(1)}.ImagePositionPatient(1);
    y_start = info{ordering(1)}.ImagePositionPatient(2);
    z_start = info{ordering(1)}.ImagePositionPatient(3);
    
    %Get modality and institution
    modality = info{ordering(1)}.Modality;
    
    switch modality
        case 'CT'
            %Get # z-slices and slice thickness
            Nz = length(list);
            CT = zeros(Nx,Ny,Nz);
            %slice thickness not in header, take different in location
            dz = info{ordering(1)}.SliceThickness;
            
        case 'MR'
            %Get # z-slices and slice thickness
            Nz = length(list);
            MR = zeros(Nx,Ny,Nz);
            %slice thickness not in header, take different in location
            dz = info{ordering(1)}.SliceThickness;
            
        case 'PT'
            %Get # z-slices and slice thickness
            Nz = double(info{ordering(1)}.NumberOfSlices);
            dz = info{ordering(1)}.SliceThickness;
            Bq = zeros(Nx,Ny,Nz);
            SUV = zeros(Nx,Ny,Nz);
    
            %SUV factor calculation
            weight = info{ordering(1)}.PatientWeight;        %in kg
            scan_datetime = ...
                [info{ordering(1)}.SeriesDate info{ordering(1)}.SeriesTime];             
            %if someone didn't enter the tracer info, then the inj_dose
            %will be blank. In that case, the script will write to
            %Bq/cc
            try
                half_life = info{ordering(1)}.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife;            %in s          
                try
                    inj_datetime = info{ordering(1)}.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartDatetime;
                catch
                    inj_date = info{ordering(1)}.StudyDate;
                    inj_time = info{ordering(1)}.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime;
                    inj_datetime = [inj_date inj_time];
                end
                inj_dose = info{ordering(1)}.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose;            %in Bq
                inj_datetimeVect = ...
                    datevec(inj_datetime, 'yyyymmddHHMMSS.FFF');
                scan_datetimeVect = ...
                    datevec(scan_datetime, 'yyyymmddHHMMSS');
                %Delay time between injection and start of scan (in s)
                t=etime(scan_datetimeVect, inj_datetimeVect);
                scan_dose = inj_dose*exp(-(log(2)*t)/half_life);%in  Bq                                            
                %Calculate SUV factor
                SUV_factor = 1/((scan_dose/weight)*0.001);
            catch err
                warning(['Tracer activity not found in DICOM' ...
                    ' header. PET image will be in units of Bq/cc']);
                inj_dose = [];
                SUV_factor = 1;
            end
    end
     
    for j = 1:Nz  
        
        if(strcmp(modality,'MR')>0)%factor to convert to correct units
        else
            m = info{ordering(j)}.RescaleSlope;
            b = info{ordering(j)}.RescaleIntercept;
        end
        
        %read data and apply scaling factors + shift
        switch modality
            case 'MR'                               
                MR(:,:,j) = ...
                    double(dicomread(info{ordering(j)}.Filename));
            case 'CT'                               
                CT(:,:,j) = ...
                    double(dicomread(info{ordering(j)}.Filename)).*m+b;
            case 'PT'         
                Bq(:,:,j) = ...
                    double(dicomread(info{ordering(j)}.Filename)).*m+b;
                SUV(:,:,j) = SUV_factor*...
                     double(dicomread(info{ordering(j)}.Filename)).*m+b;
        end
    
    end
    
    
    % OUTPUT: save to file
    Geometry = [];
    header = info{ordering(1)};
    switch modality
        case 'CT'
            Geometry.data = CT;
            Geometry.voxel_size = [dx, dy, dz];
            Geometry.start = [x_start, y_start, z_start];
            
        case 'MR'
            Geometry.data = MR;
            Geometry.voxel_size = [dx, dy, dz];
            Geometry.start = [x_start, y_start, z_start];
            
        case 'PT'
            if ~isempty(inj_dose) 
                Geometry.data = SUV;
                Geometry.SUV_factor = SUV_factor;
            else %in case someone didn't enter the dose
                Geometry.data = Bq;
            end
            Geometry.voxel_size = [dx, dy, dz];
            Geometry.start = [x_start, y_start, z_start];
                   
    end
    
end

%###############################################################################
%###############################################################################
    
    
    