function lib = calibration_fit_entry(lib, fit, tran_serial, ds_serial, isppa_wcm2, source)
% CALIBRATION_FIT_ENTRY  Library entry from the current entry and a CALIBRATION_FIT result
%
% Per-depth fits become refinements (LIBRARY_REFINEMENT_FROM_FITS). Nothing
% is written.
%
% Use as:
%   lib = calibration_fit_entry(lib, fit, tran_serial, ds_serial, isppa_wcm2, source)
%
% Input:
%   lib         - current library entry, or [] / struct() for none
%   fit         - CALIBRATION_FIT result (.model)
%   tran_serial - transducer serial for meta.tran_serial ('' = keep)
%   ds_serial   - driving-system serial for meta.driving_system ('' = none)
%   isppa_wcm2  - reference ISPPA for the stored drives [W/cm²]
%   source      - provenance string (profile CSV) for the refinements
%
% Output:
%   lib - library entry
%
% See also: CALIBRATION_FIT, CALIBRATION_SAVE_ENTRY, LIBRARY_REFINEMENT_FROM_FITS

arguments
    lib
    fit (1,1) struct
    tran_serial {mustBeTextScalar}
    ds_serial {mustBeTextScalar}
    isppa_wcm2 {mustBeNumeric}
    source {mustBeTextScalar} = ''
end

    if nargin < 6; source = ''; end
    if isempty(lib); lib = struct(); end
    if ~isempty(tran_serial); lib.meta.tran_serial = tran_serial; end
    if ~isempty(ds_serial) && ~strcmp(ds_serial, '(none)'); lib.meta.driving_system = ds_serial; end
    lib = library_refinement_from_fits(lib, fit.model, isppa_wcm2, [], source);
end
