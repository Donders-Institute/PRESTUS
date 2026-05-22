function tran_ini_data = generate_tran_ini_from_geometry(tran)
% GENERATE_TRAN_INI_FROM_GEOMETRY  Build element-position struct from ring geometry
%
% Derives the representative 3-D position of each annular element from the
% inner/outer diameters and radius of curvature stored in the transducer YAML,
% then packages them in the struct format expected by COMPUTE_PHASES.  This
% allows geometric phase computation without an external Imasonic-style .ini
% file — useful for 'generic' manufacturer transducers or as a fallback.
%
% Coordinate system (same as Imasonic .ini convention):
%   Origin  : natural focus (centre of curvature of the bowl).
%   Z-axis  : axial direction, positive toward the target / away from bowl.
%   X/Y     : lateral; y = 0 by convention for axisymmetric rings.
%   Z < 0   : all elements are behind the focal plane.
%
% For element i with ring mid-radius r_mid [mm] on a spherical bowl of
% curvature radius R [mm]:
%   x = r_mid,  y = 0,  z = -sqrt(R^2 - r_mid^2)
%
% Use as:
%   tran_ini_data = generate_tran_ini_from_geometry(tran)
%
% Input:
%   tran - transducer struct as loaded from the equipment YAML by
%          load_equipment_config(); requires fields:
%            tran.n_elem
%            tran.transducer.annular.elem_id_mm   (inner diameters [mm])
%            tran.transducer.annular.elem_od_mm   (outer diameters [mm])
%            tran.transducer.annular.curv_radius_mm
%
% Output:
%   tran_ini_data - struct with one section:
%     .elements.x1 … .xN  — position strings 'x_mm|y_mm|z_mm'
%                            identical to the format read by READ_INI_FILE
%
% See also: COMPUTE_PHASES, SET_REAL_PHASES, READ_INI_FILE

arguments
    tran (1,1) struct
end

    n_elem      = tran.n_elem;
    R_c         = tran.transducer.annular.curv_radius_mm;   % [mm]
    elem_id_mm  = tran.transducer.annular.elem_id_mm(:)';   % inner diameters [mm]
    elem_od_mm  = tran.transducer.annular.elem_od_mm(:)';   % outer diameters [mm]

    if numel(elem_id_mm) ~= n_elem || numel(elem_od_mm) ~= n_elem
        error(['generate_tran_ini_from_geometry: elem_id_mm / elem_od_mm length (%d/%d) ' ...
            'does not match n_elem (%d).'], numel(elem_id_mm), numel(elem_od_mm), n_elem);
    end

    % Mid-radius of each ring [mm]: mean of inner and outer radii
    r_mid = (elem_id_mm / 2 + elem_od_mm / 2) / 2;

    tran_ini_data.elements = struct();

    for i = 1:n_elem
        r = r_mid(i);

        if r > R_c
            error(['generate_tran_ini_from_geometry: element %d mid-radius %.2f mm exceeds ' ...
                'bowl radius of curvature %.2f mm.'], i, r, R_c);
        end

        % Z position: element sits on the bowl surface, behind the natural focus
        z = -sqrt(R_c^2 - r^2);

        field = sprintf('x%d', i);
        tran_ini_data.elements.(field) = sprintf('%.4f|0.0|%.4f', r, z);
    end

end
