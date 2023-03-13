function [sphereInds,sphereCenters]=mask2SearchlightSpheres(Mask,radius)
% Generates a cell array (sphereInds) describing searchlight spheres.

% A sphere is built around each non-zero mask entry, and hence the
% cell-array is nnz(mask) long.

% Each cell contains a list of voxels sampled by the spheres.
% Only voxels within the mask are sampled!

% The indecis within the cell are in relation to the mask non-zero entries
% (i.e. index 1 is the first non-zero element in mask).

% radius is given in voxels

% sphereCenters (optional)- x,y,z of sphere centers

% 2017 Tal Golan @ Malach Lab

mask=load_nii(Mask);
mask=mask.img;

maskInds=find(mask);
nSpheres=numel(maskInds);
[x,y,z]=ind2sub(size(mask),maskInds);
x=double(x); y=double(y); z=double(z);
sphereInds=cell(nSpheres,1);

fprintf('\nmask to spheres');
for iSphere=1:nSpheres
    
    dist=sqrt((x(iSphere)-x).^2+(y(iSphere)-y).^2+(z(iSphere)-z).^2);
    sphereInds{iSphere}=uint32(find(dist<=radius));
    
    if mod(iSphere,10000)==0
        fprintf('.');
    end
end

if nargin>1
    sphereCenters=struct;
    
    sphereCenters.x=x;
    sphereCenters.y=y;
    sphereCenters.z=z;
    
end
end