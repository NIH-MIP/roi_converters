%read in CT data
ct = ConvertDicom(['folder location of CT dicoms']);
%read in VOI data
voi = VOIreader(['name and location of VOI file']);

%initalize mask to CT size
mask = zeros(size(ct.data));
for i = 1:size(voi,1)
    slice_pts = voi{i,2}; %fill slice points based on VOI file
    slice_z = voi{i,1};
    slice_mask = poly2mask(slice_pts(:,1),slice_pts(:,2),size(t2.data,1),size(t2.data,2)); %create mask from points
    mask(:,:,slice_z) = slice_mask; %fill mask at every slice
end