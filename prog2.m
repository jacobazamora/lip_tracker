clc;
clear;
close all;

% img = imread('liptracking2_01302.jpg');
% figure;
% imshow(img);
% hold on;

load('lip_template0.mat'); % pre-made template
idx1 = 2;
idx2 = 100;

%opening the directory
srcFiles = 'liptracking2';
figure;
dinfo = dir(fullfile(srcFiles));
dinfo([dinfo.isdir]) = [];  
nfiles = length(dinfo);

%setting the active contour for the first ~100 frames
for i = idx1:idx2
 filename = fullfile(srcFiles, dinfo(i).name);
 f1 = fopen(filename, 'r');
 img = imread(filename);

 % below is tracking the mouth
 E_internal = 0;
 alpha = 0.2;
 beta = 0.003;
 lambda = 8.5;
 
 % making Matrix A
 matrix_A = zeros(length(xy), length(xy));
 
 % this for-loop fills out the matrix
 for i = 1:length(xy)
     % first row
     if (i == 1)
         matrix_A(i,1) = ((2 * alpha) + (6 * beta));
         matrix_A(i,2) = ((-alpha) - (4 * beta));
         matrix_A(i,3) = beta;
         matrix_A(i, length(xy)-1) = beta;
         matrix_A(i, length(xy)) = ((-alpha) - (4 * beta));
         % second row
     elseif (i == 2)
         matrix_A(i,1) = ((-alpha) - (4 * beta));
         matrix_A(i,2) = ((2 * alpha) + (6 * beta));
         matrix_A(i,3) = ((-alpha) - (4 * beta));
         matrix_A(i,4) = beta;
         matrix_A(i,length(xy)) = beta;
         % second-to-last row
     elseif (i == length(xy)-1)
         matrix_A(i,1) = beta;
         matrix_A(i,length(xy)-3) = beta;
         matrix_A(i,length(xy)-2) = ((-alpha) - (4 * beta));
         matrix_A(i,length(xy)-1) = ((2 * alpha) + (6 * beta));
         matrix_A(i,length(xy)) = ((-alpha) - (4 * beta));
         % last row
     elseif (i == length(xy))
         matrix_A(i,1) = ((-alpha) - (4 * beta));
         matrix_A(i,2) = beta;
         matrix_A(i,length(xy)-2) = beta;
         matrix_A(i,length(xy)-1) = ((-alpha) - (4 * beta));
         matrix_A(i,length(xy)) = ((2 * alpha) + (6 * beta));
         % for third row, third-to-last row and rows between
     else
         matrix_A(i,i-2) = beta;
         matrix_A(i,i-1) = ((-alpha) - (4 * beta));
         matrix_A(i,i) = ((2 * alpha) + (6 * beta));
         matrix_A(i,i+1) = ((-alpha) - (4 * beta));
         matrix_A(i,i+2) = beta;
     end
 end
 
 % making the identity matrix
 I_matrix = eye(length(xy));
 
 % finding the image force
 % we will need a grayscale image to find the gradients
 grayImg = rgb2gray(img);

 gausImg = imgaussfilt(grayImg);
 doubleImg = double(gausImg);
 [Gx,Gy] = gradient(doubleImg);% [gradient of image in respect to x, graident of image in respect to y]
 [Gx_x,Gx_y] = gradient(Gx); % [gradient of Ix in respect to x, gradient of Ix in respect to y]
 [Gy_x,Gy_y] = gradient(Gy); % [gradient of Iy in respect to x, gradient of Iy in respect to y]
 
 f_x_xy = zeros(length(xy),1); % Ax + f_x_(x,y) = 0
 f_y_xy = zeros(length(xy),1); % Ax + f_y_(x,y) = 0
 for i = 1:length(xy)
     f_x_mat = 2 * ((Gx(floor(xy(2,i)),floor(xy(1,i))) * Gx_x) + (Gy(floor(xy(2,i)),floor(xy(1,i))) * Gy_x));
     f_y_mat = 2 * ((Gx(floor(xy(2,i)),floor(xy(1,i))) * Gx_y) + (Gy(floor(xy(2,i)),floor(xy(1,i))) * Gy_y));
     f_x_xy(i,1) = f_x_mat(floor(xy(2,i)),floor(xy(1,i)));
     f_y_xy(i,1) = f_y_mat(floor(xy(2,i)),floor(xy(1,i)));
 end
 
 % for j = 1:5
 new_x = zeros(length(xy),1); % new x values after every iteration
 new_y = zeros(length(xy),1); % new y values after every iteration
 % for i = 1:length(xy(1,:))
 %     new_x(i,1) = inv(matrix_A+(lambda*I_matrix))*((lambda * transpose(xy(1,i)))-f_x_xy(i,1));
 %     new_y(i,1) = inv(matrix_A+(lambda*I_matrix))*((lambda * transpose(xy(2,i)))-f_y_xy(i,1));
 % end
 xy_t=transpose(xy);
 new_x=inv(matrix_A+(lambda*I_matrix))*(lambda*xy_t(:,1)-f_x_xy);
 new_y=inv(matrix_A+(lambda*I_matrix))*(lambda*xy_t(:,2)-f_y_xy);

 for i = 1:length(xy)
     xy(1,i) = new_x(i,1);
     xy(2,i) = new_y(i,1);
 end
 % end
 
 % plotting new contour points
 imshow(img);
 hold on;
 xys = spline(t,xy,ts);
 xs = xys(1,:); 
 ys = xys(2,:);
 plot(xs, ys, 'b-',xy(1,:),xy(2,:),'r*');
 drawnow;
end
% fclose(f1);
