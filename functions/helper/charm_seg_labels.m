function labels = charm_seg_labels()
% CHARM_SEG_LABELS  Return the fixed SimNIBS charm tissue label indices
%
%   Returns a struct with named fields mapping tissue names to their integer
%   label values as produced by SimNIBS charm segmentation. These values
%   are fixed by the charm pipeline and should not need to be changed.
%
% Use as:
%   labels = charm_seg_labels()
%
% Output:
%   labels - [1x1] struct with fields:
%     bonemask        - All tissues within the skull [1,2,3,4,7,8,9]
%     intracranial    - Intracranial tissues (skull excluded) [1,2,3,9]
%     external        - Background / external [0]
%     wm              - White matter [1]
%     gm              - Grey matter [2]
%     csf             - Cerebrospinal fluid [3]
%     skull           - Skull (combined) [4]
%     skin            - Skin [5]
%     eye             - Eye [6]
%     skull_cortical  - Cortical bone [7]
%     skull_trabecular- Trabecular bone [8]
%     blood           - Blood vessels [9]
%     muscle          - Muscle [10]
%
% See also: GETIDX, CHECK_LAYERS

labels.bonemask         = [1, 2, 3, 4, 7, 8, 9];
labels.intracranial     = [1, 2, 3, 9];
labels.external         = 0;
labels.wm               = 1;
labels.gm               = 2;
labels.csf              = 3;
labels.skull            = 4;
labels.skin             = 5;
labels.eye              = 6;
labels.skull_cortical   = 7;
labels.skull_trabecular = 8;
labels.blood            = 9;
labels.muscle           = 10;

end
