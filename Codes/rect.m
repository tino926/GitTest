clear all;
tic
% l = [253, 132, 259, 210; 373 18 371 126; 541 50 522 197];
% img = double(imread('.\2_640x480\C1_000000328.bmp'))/255;
% l = [320, 373, 338, 481; 34 560 80 665; 1149 279 1079 498];
% img = double(imread('pic.JPG'))/255;
% l = [111, 95, 172, 412; 383 150 379 274; 482 189 467 306];
% img = double(imread('Picture 15.JPG'))/255;
l = [227, 133, 251, 358; 513 134 497 346];
img = double(imread('syn01.png'))/255;
h = size(img,1);
w = size(img,2);

for i = 1:size(l,1);
    % x = ay+b;
    a(i) = (l(i, 3) - l(i,1)) / (l(i,4) - l(i,2));
    b(i) = l(i,1) - l(i,2) * a(i);
end

% x - y*a1 = b1 ...
% [1 -a1] [x, y]' = b1
X = [ ones(length(a),1) -a'] \ b';

tmp = ( [(1:w)' ones(w,1)  ; (1:w)' ones(w,1)*h ; ...
    ones(h,1) (1:h)'; ones(h,1)*w (1:h)'] - repmat(X',w+w+h+h,1));
dist2X = sqrt(sum(tmp.^2,2));
tmp = ([1 1; w 1; 1 h; w h] - repmat(X',4,1));
ang2X = -atan(tmp(:, 1)./tmp(:,2));

maxDist = max(dist2X);
minDist = min(dist2X);
maxAng = max(ang2X);
minAng = min(ang2X);

% scale(y) = (maxDist-y+1)/maxDist
% S(y) = maxDist/(maxDist-y+1);
% integral of S(y) =  maxDist * (-ln(maxDist-y+1))
S_0 = maxDist * (-log(maxDist-0+1));
h_n_ = ...
    ceil(maxDist * (-log(maxDist-  (maxDist-minDist+1)  +1)) - S_0);


h_n = ceil(maxDist-minDist+1);
w_n = ceil((maxAng - minAng)*maxDist);

img_rect = zeros(h_n, w_n, 3);
img_rect_ = zeros(h_n_, w_n, 3);

x2Ang = zeros(1,w_n);
for x = 1:w_n
    x2Ang(x) = (x-1)/(maxDist-1) + minAng;
end
y2Dist = zeros(1,h_n);
for y = 1:h_n
    y2Dist(y) = maxDist-y+1;
end
y2Dist_ = zeros(1,h_n_);
for y = 1:h_n_
    y2Dist_(y) = maxDist+1 - ...
        (-(exp(-(y+S_0)/maxDist)-maxDist-1));
end
toc

tic
for x = 1:w_n
    p_Ang = x2Ang(x);
    for y = 1:h_n
        p_dist = y2Dist(y);
        
        x_o = p_dist*sin(p_Ang) + X(1);
        y_o = X(2) - p_dist*cos(p_Ang);
        
        x_o_p = floor(x_o);
        x_o_n = ceil(x_o);
        y_o_p = floor(y_o);
        y_o_n = ceil(y_o);
        
        if (x_o_p > 0 && x_o_n <=w && y_o_p >0 && y_o_n <=h)
%             img_rect(y,x,:) = img(y_o,x_o,:);
            dx_o_p = x_o - x_o_p;
            dx_o_n = x_o_n - x_o;
            dy_o_p = y_o - y_o_p;
            dy_o_n = y_o_n - y_o;
            
            img_rect(y,x,:) = ...
                img(y_o_p, x_o_p, :) * dy_o_n * dx_o_n+ ...
                img(y_o_p, x_o_n, :) * dy_o_n * dx_o_p + ...
                img(y_o_n, x_o_p, :) * dy_o_p * dx_o_n + ...
                img(y_o_n, x_o_n, :) * dy_o_p * dx_o_p;
        end
    end
    
    for y = 1:h_n_
        p_dist = y2Dist_(y);
        
        x_o = p_dist*sin(p_Ang) + X(1);
        y_o = X(2) - p_dist*cos(p_Ang);
        
        x_o_p = floor(x_o);
        x_o_n = ceil(x_o);
        y_o_p = floor(y_o);
        y_o_n = ceil(y_o);
        
        if (x_o_p > 0 && x_o_n <=w && y_o_p >0 && y_o_n <=h)
%             img_rect(y,x,:) = img(y_o,x_o,:);
            dx_o_p = x_o - x_o_p;
            dx_o_n = x_o_n - x_o;
            dy_o_p = y_o - y_o_p;
            dy_o_n = y_o_n - y_o;
            
            img_rect_(y,x,:) = ...
                img(y_o_p, x_o_p, :) * dy_o_n * dx_o_n+ ...
                img(y_o_p, x_o_n, :) * dy_o_n * dx_o_p + ...
                img(y_o_n, x_o_p, :) * dy_o_p * dx_o_n + ...
                img(y_o_n, x_o_n, :) * dy_o_p * dx_o_p;
        end
    end
end
toc


for i = 10:30:size(img,2)-1
    a_tmp = (i - X(1)) / (1 - X(2));
    b_tmp = i - a_tmp;
    
    for j = 1:size(img,1);
        img( j, round(a_tmp*j + b_tmp), 1) = 255;
    end
end

