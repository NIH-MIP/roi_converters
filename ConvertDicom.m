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
