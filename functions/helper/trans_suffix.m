function s = trans_suffix(ti, n_trans)
% TRANS_SUFFIX  Return '_trans-XX' when n_trans > 1, otherwise ''.
%
%   s = trans_suffix(ti, n_trans)
%
%   ti      - transducer index (1-based)
%   n_trans - total number of transducers in the run
%
%   When n_trans == 1 the suffix is omitted so single-transducer outputs
%   stay uncluttered.  When n_trans > 1 the suffix '_trans-XX' is appended
%   (zero-padded to two digits) to disambiguate per-transducer files.
%
% Example
%   trans_suffix(2, 3)  % '_trans-02'
%   trans_suffix(1, 1)  % ''

if n_trans > 1
    s = sprintf('_trans-%02d', ti);
else
    s = '';
end
end
