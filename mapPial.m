function pial=mapPial(region)

VG=spm_vol(region(1,:));
pial=zeros(VG.dim(1:3)); 
for i=1:VG.dim(3),
  pial(:,:,i) = spm_slice_vol(VG,spm_matrix([0 0 i]),VG.dim(1:2),1);
end

end