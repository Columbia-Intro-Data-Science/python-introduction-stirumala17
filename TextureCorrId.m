%Sragvi Tirumala
%Texture test

%Texture Research
close all
clc;

%glcm vars
offsets = [0 1; -1 1;-1 0;-1 -1];
numoffsets = length(offsets);

%Reading images from folder and storing it in variable
path = uigetdir('Final_Images/Final_Images');
imdir = dir(path);

%store glcm stats
diseasestats = cell(1,length(imdir)-2);
diseasenames = cell(1,7);

%number of images used to gauge data
numimages = 20;

for i = 1:length(imdir)-2
    disease = imdir(i+2).name;
    diseasenames{i} = disease;
    dispath = [path,'/',disease];
    disdir = dir(dispath);
    
    
    diseasestats{i} = cell(1,5);
    for j = 1:numimages
        subject = disdir(j+2).name;
        img = imread([dispath,'/',subject]);
        img = rgb2gray(img);
        
        glcm = graycomatrix(img,'Offset',offsets,'Symmetric',true);
        
        stats = graycoprops(glcm);
        diseasestats{i}{j} = stats;
    end
end

%Plotting all the goddamn stats
contrast = zeros(numimages,numoffsets,7);
correlation = zeros(numimages,numoffsets,7);
energy = zeros(numimages,numoffsets,7);
homogeneity = zeros(numimages,numoffsets,7);

for i = 1:length(imdir)-2
    for j = 1:numimages
        contrast(j,:,i) = diseasestats{i}{j}.Contrast;
        correlation(j,:,i) = diseasestats{i}{j}.Correlation;
        energy(j,:,i) = diseasestats{i}{j}.Energy;
        homogeneity(j,:,i) = diseasestats{i}{j}.Homogeneity;
    end
end

%average stats over images used for each disease
MCon = mean(contrast,1);
MCor = mean(correlation,1);
ME = mean(energy,1);
MH = mean(homogeneity,1);

SCon = std(contrast,1);
SCor = std(correlation,1);
SE = std(energy,1);
SH = std(homogeneity,1);

c = categorical({'ato','pit','con','mel','pso','rin','mim'});

% Contrast

figure
subplot(2,2,1)
errorbar(c,reshape(MCon(1,1,:),[1,7]),reshape(SCon(1,1,:),[1,7]))
title('Offset 0,1')

subplot(2,2,2)
errorbar(c,reshape(MCon(1,2,:),[1,7]),reshape(SCon(1,2,:),[1,7]))
title('Offset -1,1')

subplot(2,2,3)
errorbar(c,reshape(MCon(1,3,:),[1,7]),reshape(SCon(1,3,:),[1,7]))
title('Offset -1,0')

subplot(2,2,4)
errorbar(c,reshape(MCon(1,4,:),[1,7]),reshape(SCon(1,4,:),[1,7]))
title('Offset -1,-1')

suptitle('Contrast')

% Correlation

figure
subplot(2,2,1)
bar(c,reshape(MCor(1,1,:),[1,7]))
title('Offset 0,1')

subplot(2,2,2)
bar(c,reshape(MCor(1,2,:),[1,7]))
title('Offset -1,1')

subplot(2,2,3)
bar(c,reshape(MCor(1,3,:),[1,7]))
title('Offset -1,0')

subplot(2,2,4)
bar(c,reshape(MCor(1,4,:),[1,7]))
title('Offset -1,-1')

suptitle('Correlation')

%Energy

figure
subplot(2,2,1)
bar(c,reshape(ME(1,1,:),[1,7]))
title('Offset 0,1')

subplot(2,2,2)
bar(c,reshape(ME(1,2,:),[1,7]))
title('Offset -1,1')

subplot(2,2,3)
bar(c,reshape(ME(1,3,:),[1,7]))
title('Offset -1,0')

subplot(2,2,4)
bar(c,reshape(ME(1,4,:),[1,7]))
title('Offset -1,-1')

suptitle('Energy')

% Homogeneity

figure
subplot(2,2,1)
bar(c,reshape(MH(1,1,:),[1,7]))
title('Offset 0,1')

subplot(2,2,2)
bar(c,reshape(MH(1,2,:),[1,7]))
title('Offset -1,1')

subplot(2,2,3)
bar(c,reshape(MH(1,3,:),[1,7]))
title('Offset -1,0')

subplot(2,2,4)
bar(c,reshape(MH(1,4,:),[1,7]))
title('Offset -1,-1')

suptitle('Homogeneity')







