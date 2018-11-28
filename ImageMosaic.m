%% Use SIFT to get matches between images
img1 = imread('test2_1.jpg');
img2 = imread('test2_2.jpg');

[f1, d1] = vl_sift(single(rgb2gray(img1))); % [x, y, s, th]
[f2, d2] = vl_sift(single(rgb2gray(img2)));
matches = vl_ubcmatch(d1,d2); % [index in f1, index in f2]


%% Use RANSAC for homography computation
% Iterate N times
N = 100; % number of times to iterate
t = 3; % threshold for whether pixel matches
count = 0; % number of matches for currH
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
    
    if numInliers>count
        count = numInliers;
        currH = H;
    end
end