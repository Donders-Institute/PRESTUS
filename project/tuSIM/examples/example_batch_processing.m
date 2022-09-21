addpath '/home/common/matlab/fieldtrip/qsub'

% example batch calibration
source_amp_arr = 60000:10000:140000;
source_phase_deg = [0.0 71.6 143.2 214.8; 0.0 355.3 350.6 341.1];

source_amp_arr = 93000:1000:97000;
source_phase_deg = [0.0 355.3 350.6 341.1];

jobs = {};
for j = 1:size(source_phase_deg, 1)
    for i = 1:length(source_amp_arr)
        jobs{length(jobs)+1} = qsubfeval(@calibrate_in_water, 'sjoerd_config.yaml', source_amp_arr(i), source_phase_deg(j,:), 'timreq',  60*10,  'memreq',  20*(1024^3),  'options', '-l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=5.0"');
    end
end

% example runnning several subjects for left and right amygdala with
% 'pipeline_sjoerd'

subject_ids = [1,5,6];

jobs = {};

for subject_id = [1,5,6]
    for target_id = [1,2]
        jobs{length(jobs)+1} = qsubfeval(@pipeline_sjoerd, subject_id, target_id, 'timreq',  60*60*7,  'memreq',  50*(1024^3),  'options', '-l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=5.0"');
    end
end


% example runnning several subjects for with 'example_pipeline_judith'
% single session

% batch jobs
for subject_id = [1, 12]
    qsubfeval(@example_pipeline_judith, subject_id, 'timreq',  60*60*7,  'memreq',  50*(1024^3),  'options', '-l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=8.0"');
end
