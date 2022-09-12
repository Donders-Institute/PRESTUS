function single_subject_pipeline_wrapper(subject_id, parameters)

% single_subject_pipeline_wrapper adds paths to functions and toolboxes and runs single_subject_pipeline

    addpath('functions')
    addpath(genpath('toolboxes')) % add toolboxes

    single_subject_pipeline(subject_id, parameters);

end

