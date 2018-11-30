%% Use SIFT to get matches between images
img1 = imread('test1_1.png');
img2 = imread('test1_2.png');

[f1, d1] = vl_sift(single(rgb2gray(img1))); % f: [x, y, s, th]
[f2, d2] = vl_sift(single(rgb2gray(img2)));
matches = vl_ubcmatch(d1,d2); % [index in f1, index in f2]
numMatches = size(matches,2); % number of matches

%% Use RANSAC for homography computation
% Iterate N times
N = 100; % number of times to iterate
t = 3; % threshold for whether pixel matches
allCount = {}; % array to hold all inliers
allH = {}; % array to hold all Hs
count = 0; % number of matches for currH
maxInliers = 0; % maximum inliers
currH = zeros(3,3); % H with most inliers - current best
for n=1:N
    % Randomly select a sample
    index=randi(size(matches,2));
    match=matches(:,index);
    
    % Compute homography with Direct Linear Transformation
    x1 = f1(1, match(1));
    y1 = f1(2, match(1));
    x2 = f2(1, match(2));
    y2 = f2(2, match(2));
    deltaX = x1-x2;
    deltaY = y1-y2;
    H=[1 0 deltaX; 0 1 deltaY; 0 0 1];
    
    % Project points from x to x' for each potential match
    numInliers = 0;
    for x_ind = 1:size(matches,2)
        % x' = H*x
        newMatch = matches(:,x_ind);
        x_f1 = f1(1,newMatch(1));
        y_f1 = f1(2,newMatch(1));
        x_prime = H * [x_f1; y_f1; 1];  % from img1
        
        % Calculate distance/error
        x_f2 = f2(1,newMatch(2));
        y_f2 = f2(2,newMatch(2));
        error = (x_prime(1)-x_f2)^2 + (x_prime(2)-y_f2)^2;
        
        % Check if under threshold
        if error<t
            numInliers = numInliers+1;
        end
    end
    
    % Save current homography's inliers
    allCount{n} = numInliers;
    allH{n} = H;
    
    
    % Update maximum number of inliers
    if numInliers>count
        count = numInliers;
        currH = H;
    end
end

% % Get the best H and number of inliers
% [allCount,ind] = max(allCount);
% H = H{ind};
% maxInliers = allCount{ind};
% sprintf('%d max inliers', maxInliers);

%% Visualize all matches
% Pad if widths of two images are different
width_diff1 = max(size(img2,1)-size(img1,1),0);
width_diff2 = max(size(img1,1)-size(img2,1),0);

% Plot match lines
figure();
imagesc([padarray(img1,width_diff1,'post') padarray(img2,width_diff2,'post')]);
img1_height = size(img1,2);
line([f1(1,matches(1,:));f2(1,matches(2,:))+img1_height],[f1(2,matches(1,:));f2(2,matches(2,:))]);
title(sprintf('%d matches', numMatches));
axis image off;

drawnow;

%% Visualize inlier matches
figure();
imagesc([padarray(img1,width_diff1,'post') padarray(img2,width_diff2,'post')]);
img1_height = size(img1,2);
line([f1(1,matches(1,count));f2(1,matches(2,count))+img1_height],[f1(2,matches(1,count));f2(2,matches(2,count))]);
title(sprintf('%d inliner matches', count));
axis image off;

drawnow ;