function [] = voiMaker(saveLoc, structures)   
    for k=1:length(structures)
        roi = structures{k,1}; %grab single ROI and write to seperate folder
        load('voi_mainBase.mat'); % i already copied a base example to use as the structure
        voi = voi_mainBase;
        fileID = fopen([saveLoc roi.name '.txt'],'w'); %write as text file, then at the end we copy and make to VOI
        
        %write base
        for mm = 1:size(voi_mainBase,1)
            fprintf(fileID,'%s\t%s\t%s\r\n',voi_mainBase{mm,1},voi_mainBase{mm,2},voi_mainBase{mm,3});
        end
        
        %write points for ROI based on slice number
        %this iterates over the number of SLICES which are saved as
        %seperate arrays within rtstruct file
        for i = 1:size(roi.points,2)
            load('voi_sliceBase.mat');
            sliceDat = roi.points{1,i};
            slice = unique(sliceDat(:,3));
            xy = sliceDat(:,1:2);

            voi_sliceBase{1,1} = int2str(slice-1); %VOIs subtract one for slice number
            voi_sliceBase{2,1} = int2str(1); %dummy variable
            voi_sliceBase{3,1} = int2str(size(xy,1)); %number of points in VOI

            %copy over slice base for each slice and fill w variables
            for mm = 1:size(voi_sliceBase,1)
                fprintf(fileID,'%s\t%s\t%s\r\n',voi_sliceBase{mm,1},voi_sliceBase{mm,2},voi_sliceBase{mm,3});
            end

            %then 
            for j = 1:size(xy,1)
                voi_sliceBase{j+3,1} = [xy(j,1) xy(j,2)];
                fprintf(fileID,'%0.3f %0.3f\r\n',xy(j,1),xy(j,2));
            end

            voi = cat(1,voi,voi_sliceBase);
            clear voi_sliceBase
        end

        voi_unique = int2str(k);
        voi_unique_text = '# unique ID of the VOI';
        fprintf(fileID,'%s\t%s\t%s\r\n',voi_unique,'',voi_unique_text);
        
        fclose(fileID);
        copyfile([saveLoc roi.name '.txt'],[saveLoc roi.name '.voi']);
        delete([saveLoc roi.name '.txt']);
    end
end