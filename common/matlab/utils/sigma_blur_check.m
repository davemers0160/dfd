format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);

%% start of main code
maxSigma = 3;
classes = 256;
point = zeros(50,50);
kernel_size = 43;

point(24:26, 24:26) = 255;

BlurStep = maxSigma / classes;

sigma=[];
c=[];

for idx=1:classes
%     sigma(idx) = 2*(((idx)/classes)^0.3);
    sigma(idx) = BlurStep*idx;

    %sigma(idx) = (maxSigma*2*BlurStep*idx)/(1 + 2 * abs(BlurStep*idx));
    %sigma(idx) = (maxSigma*0.2*BlurStep*idx)/(maxSigma-0.8*abs(BlurStep*idx));
    
    kernel = createGaussKernel(kernel_size, sigma(idx));

    b_tmp = imfilter(point, kernel);
    %b_tmp = imgaussfilt(point, sigma(idx));
    
    B = floor(b_tmp);
    
    b2 = B(1:floor(size(point,1)/2),floor(size(point,1)/2));
    c(idx) = sum(b2>0)-1;
end


%% plot results

figure(1);
plot(1:classes, sigma,'.-b');



figure(2);
plot(1:classes, c, '-g');