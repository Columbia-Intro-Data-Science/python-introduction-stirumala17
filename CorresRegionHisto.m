%% This HistogramCompare program is used for finding abnormal temperature region on the corresponding side. 
%  GWU BME MedImg Lab --- Advisor : Prof. Murray Loew
%  Author : Shijian Fan
%  Verson_1

%% Pre- Processing
close all
clear all
clc

%Read the original image

ImgFile = '0040DN.tif';
OriginalImg  = imread(ImgFile);

%figure, imshow(OriginalImg);

% Denoise using median filter

BreastImg = medfilt2(OriginalImg,[3 3],'symmetric');
figure,imshow(histeq(BreastImg));


% Pick up the reference point
[ReferPointX,ReferPointY]=ginput(2);
%Store the Reference Position
ReferLeft = [ReferPointX(1),ReferPointY(1)];
ReferRight = [ReferPointX(2),ReferPointY(2)];


%% First we find the local warm region on either side of the Breast.

%Here we selected the test region and add some noise on it

Rect = [290 142 74 74] ;
SelectedImage = imcrop(BreastImg,Rect);
%Adding white noise to the selected region(normalized the parameters)

%figure, imshow(SelectedImage);
% For imnoise function, you must use normalized noise parameters. 
OneMask = mat2gray(zeros(size(SelectedImage)));
OneMask = im2uint16(OneMask);
AddingNoise = imnoise(OneMask,'gaussian',((mean(SelectedImage(:))+400)/65535),0.001);
%figure,imshow(AddingNoise);
TestImage = BreastImg;
for i = Rect(2) : (Rect(2)+Rect(4))
    for j = Rect(1) : (Rect(1)+Rect(3))
        TestImage(i,j) = AddingNoise(i-(Rect(2)-1),j-(Rect(1)-1));
    end
end
figure, imshow(histeq(TestImage));
hold on;


%% Find the large corresponding region on the other side
%Plot a line between to reference points
plot([ReferLeft(1),ReferRight(1)],[ReferLeft(2),ReferRight(2)],'--gs','Color','b','LineWidth',0.75,'MarkerSize',10);

% The angel between two reference points
U = [1 0];
Vector_between_ReferencePoints = [ (ReferRight(1)-ReferLeft(1)) (ReferRight(2)-ReferLeft(2)) ];
Cos_Ang_Shift = dot(U,Vector_between_ReferencePoints)/(norm(U)*norm(Vector_between_ReferencePoints));
Ang_Shift = acos(Cos_Ang_Shift)*180/pi;

% The parameters that help find the corresponding region

TargetLeftX = Rect(1)+(Rect(3)/2);
TargetLeftY = Rect(2)+(Rect(4)/2);
Vector_between_TargetLeft = [ (TargetLeftX-ReferLeft(1)) (TargetLeftY-ReferLeft(2)) ];
Cos_Delta_Refer = dot(U,Vector_between_TargetLeft)/(norm(U)*norm(Vector_between_TargetLeft));
Delta_Refer = acos(Cos_Delta_Refer)*180/pi;
Delta = Delta_Refer - Ang_Shift;
Distance = sqrt((TargetLeftX - ReferLeft(1)).^2 + (TargetLeftY - ReferLeft(2)).^2);

% Plot the Target point and the reference line

plot([ReferLeft(1),TargetLeftX],[ReferLeft(2),TargetLeftY],'-gs','Color','b','LineWidth',0.75,'MarkerSize',10);
%plot(TargetLeftX,TargetLeftY,'x','LineWidth',2,'Color','b');

% Calculate the position of the right corresponding point

TargetRightX = ReferRight(1)-cosd(Delta-Ang_Shift)*Distance;
TargetRightY = ReferRight(2)-sind(Delta-Ang_Shift)*Distance;

% Plot the Target point and the reference line

plot([ReferRight(1),TargetRightX],[ReferRight(2),TargetRightY],'-gs','Color','b','LineWidth',0.75,'MarkerSize',10);
%plot(TargetRightX,TargetRightY,'x','LineWidth',2,'Color','yellow');

% Define the width and height of the corresponding region

Width = Rect(3)+0.25*Rect(3);
Height = Rect(4)+0.25*Rect(4);
Rect2 =[(TargetRightX-(Width/2)) (TargetRightY-(Height/2)) Width Height];

Corresponding_Region = imcrop(BreastImg,Rect2);

% Plot the larger Corresponding region

plot([(TargetRightX-(Width/2)) (TargetRightX-(Width/2))],[(TargetRightY-(Height/2)) (TargetRightY+(Height/2))],'-gs','Color','b','LineWidth',0.75,'MarkerSize',10);
plot([(TargetRightX-(Width/2)) (TargetRightX+(Width/2))],[(TargetRightY-(Height/2)) (TargetRightY-(Height/2))],'-gs','Color','b','LineWidth',0.75,'MarkerSize',10);
plot([(TargetRightX-(Width/2)) (TargetRightX+(Width/2))],[(TargetRightY+(Height/2)) (TargetRightY+(Height/2))],'-gs','Color','b','LineWidth',0.75,'MarkerSize',10);
plot([(TargetRightX+(Width/2)) (TargetRightX+(Width/2))],[(TargetRightY-(Height/2)) (TargetRightY+(Height/2))],'-gs','Color','b','LineWidth',0.75,'MarkerSize',10);


figure, imshow(Corresponding_Region);


%%Calculate the histogram of the selected region

 Bins =double(min(Corresponding_Region(:)))-10 : 1 : double(max(Corresponding_Region(:)))+10;
HistSelected = histc(reshape(Corresponding_Region,[1,size(Corresponding_Region,1)*size(Corresponding_Region,2)]), Bins);
figure,bar(Bins,HistSelected,'b');
ylim([0 100]);


%% Find the pixels within the A% highest value region of the histogram.
A = 1;
NumberOfPixels = size(Corresponding_Region,1)*size(Corresponding_Region,2); 
Count = 0; 
Position = size(HistSelected,2);
for i = size(HistSelected,2) : -1 : 1
    Count = Count + HistSelected(i);
    if Count >= (A*0.01*NumberOfPixels)
        Position = i;
        break
    end
end

Distance = size(HistSelected,2)- Position;
Threshold = max(Corresponding_Region(:))-Distance;
[PositionY,PositionX] = ind2sub(size(Corresponding_Region),find(Corresponding_Region>=Threshold));

DetecedImage = zeros(size(Corresponding_Region));
% Highlight the abnormal tempeture region
DetecedImage(PositionY,PositionX) = Corresponding_Region(PositionY,PositionX) ;
DetecedImage =im2uint16(mat2gray(DetecedImage));
figure, imshow(DetecedImage,[]);
%plot(PositionX,PositionY ,'r','Markersize','2'); plot lines instead

%% Calculate the Maximum/Minimum distance between those pixels 

%Merge the two distance into one matrix
ObservationMatrix = [PositionX PositionY];
DistanceBetweenPixels = pdist(ObservationMatrix);
MaxDistance = max(DistanceBetweenPixels)






