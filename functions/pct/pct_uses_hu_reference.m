function tf = pct_uses_hu_reference(parameters)
% PCT_USES_HU_REFERENCE  Whether a pCT mapping of this run scales by the densest-bone HU
%
% True when pct.mapping_density is 'marsac-babelbrain' or 'aubry',
% pct.mapping_soundspeed is 'aubry', or pct.mapping_attenuation is 'mueller'
% or 'aubry' (PCT_HU_REFERENCE), or when the BabelBrain backend maps HU itself
% (babelbrain.ct_mapping).
%
% Use as:
%   tf = pct_uses_hu_reference(parameters)
%
% See also: PCT_HU_REFERENCE, MEDIUM_PCT_HU_MAX

    arguments
        parameters (1,1) struct
    end
    m = @(key) char(get_field(parameters, {'pct', key}, 'none'));
    tf = ismember(m('mapping_density'), {'marsac-babelbrain', 'aubry'}) || ...
         strcmp(m('mapping_soundspeed'), 'aubry') || ...
         ismember(m('mapping_attenuation'), {'mueller', 'aubry'}) || ...
         (strcmp(char(get_field(parameters, {'acoustic_backend'}, '')), 'babelbrain') && ...
          ~isempty(char(get_field(parameters, {'babelbrain', 'ct_mapping'}, ''))));
end
