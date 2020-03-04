%
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%
%% Script for biplane registration
%  input: qd1 and qd2
%  output: tform (affine parameters)

%%
function tform = biplane_registration(qd1,qd2)


fixed = max(qd1,[],3);
moving = max(qd2,[],3);


fixed = fixed ./ max(max(fixed));
moving = moving ./ max(max(moving));

rate = mean(mean(fixed)) / mean(mean(moving));
moving = moving * rate;


optimizer = registration.optimizer.RegularStepGradientDescent;
metric = registration.metric.MeanSquares;
%affine
tform = imregtform(moving, fixed, 'affine', optimizer, metric);

if tform.T(1,1) < 0.9 || tform.T(2,2) < 0.9  % if scale too much, didn't registrate
    tform.T = [1 0 0; 0 1 0; 0 0 1];
end

% movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
% h = joinchannels('RGB', fixed, movingRegistered);
% dipshow(h)
% figure; imshowpair(fixed, movingRegistered,'ColorChannels',[1 2 0], 'Scaling','joint');
% axis tight
% title('Registrated image');


