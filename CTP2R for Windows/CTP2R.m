clear all

%****************Read in data**************************************
[FileName2,PathName2] = uigetfile('*.tif','Select the Image file');
Folder_Name = uigetdir(PathName2,'Select the TracksMeasurement CSV folder');

prompt = {'px -> micron conversion (xy):','px -> micron conversion (z):','Start of timelapse (t0)','Timelapse interval (delta t)'};
dlg_title = 'Image properties';
num_lines = 1;
defaultans = {'0.62','1','5','0.5'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

if ismac == 1
    sep = '/';
else if ispc == 1
        sep = '\';
    end
end

num1 = importfile([Folder_Name sep 'sheet0.csv']);
        
no_tracks = length(num1(:,2)); %this is the total number of individual "cells" that have been tracked
no_frames = max(num1(:,3)) + 1;
num_all = zeros(no_tracks,no_frames,50);

filelist = dir(Folder_Name);

for i = 1:no_tracks;
    
    num = importfile1([Folder_Name sep 'sheet' num2str(i) '.csv']);

    num_all(i,1:size(num,1),1:size(num,2)) = num;
    
end

num_all(:,:,[1:2,7:10,14:15,28:34,37:end]) = [];


num_all(find(num_all==0)) = NaN;
for i = 1:no_tracks;
if isnan(num_all(i,1,4))
    num_all(i,1,4) = 0;
end
end

%convert xy to microns
conversion = str2num(answer{1}); %px/um
num_all(:,:,[1:2,5:6]) = num_all(:,:,[1:2,5:6])*conversion; %only pos and size in px, rest in um
num_all(:,:,[3,7]) = num_all(:,:,[3,7])*str2num(answer{2}); %convert z to um

time = linspace(str2num(answer{3}),str2num(answer{4})*(no_frames-1) + str2num(answer{3}),no_frames);
time_mat = repmat(time,no_tracks,1);
deltat = linspace(0,str2num(answer{4})*(no_frames-2),no_frames-1);
deltat_mat = repmat(deltat,no_tracks,1);
deltat2 = linspace(0,str2num(answer{4})*(no_frames-3),no_frames-2);
deltat2_mat = repmat(deltat2,no_tracks,1);

%*************Timescales analysis*************************
%MSD calculation
msd_all = zeros(no_tracks,no_frames-1);

for cell = 1:no_tracks;
%calculate MSD for each trajectory
    for d_t = 1:no_frames-1;
     msd_sum = 0;
        flag = 0;
        for i = 1:no_frames-d_t;
            if isnan(num_all(cell,i+d_t,1));
                break
            else
                msd = (num_all(cell,i,1)-num_all(cell,i+d_t,1))^2 + (num_all(cell,i,2)-num_all(cell,i+d_t,2))^2 + (num_all(cell,i,3)-num_all(cell,i+d_t,3))^2;
                msd_sum = msd_sum + msd;
                flag = flag + 1;
            end
        end
     msd_d_t(d_t) = msd_sum/flag;
    end
    msd_all(cell,:) = msd_d_t;
end
      
%Directional autocorrelation calculation
da_all = zeros(no_tracks,no_frames-2);

for cell = 1:no_tracks;
%calculate autocorrelation for each trajectory
    for d_t = 1:no_frames-2;
     da_sum = 0;
        flag = 0;
        for i = 1:no_frames-1-d_t;
            if isnan(num_all(cell,i+d_t+1,1));
                break
            else
                %find both vectors and normalise
                v1 = [num_all(cell,i+1,1) - num_all(cell,i,1),num_all(cell,i+1,2)-num_all(cell,i,2),num_all(cell,i+1,3)-num_all(cell,i,3)];
                v2 = [num_all(cell,i+1+d_t,1) - num_all(cell,i+d_t,1),num_all(cell,i+1+d_t,2)-num_all(cell,i+d_t,2),num_all(cell,i+1+d_t,3)-num_all(cell,i+d_t,3)];
                v1 = v1/norm(v1);
                v2 = v2/norm(v2);
                da = dot(v1,v2);
                da_sum = da_sum + da;
                flag = flag + 1;
            end
        end
     da_d_t(d_t) = da_sum/flag;
    end
    da_all(cell,:) = da_d_t;
end

%************Split off myotome vs. non-myotome******************
% checking track positions for split tracks and position in myotome
fname = [PathName2 FileName2];
info = imfinfo(fname);
no_images = numel(info); 
im = zeros(info(1).Height, info(1).Width, no_images/no_frames);
for k = 1:no_images/no_frames
    im(:,:,k) = imread(fname, k, 'Info', info);
end

%Manual selection of ROIs
prompt = {'Select number of ROIs'};
dlg_title = 'ROI selection';
num_lines = 1;
defaultans = {'1'};
roi_ans = inputdlg(prompt,dlg_title,num_lines,defaultans);
no_roi = str2num(roi_ans{1});

max_proj = max(im,[],3);

strflag = 'No'; %repeats loop until dialog is selected as okay
while strcmp(strflag,'Yes') ~= 1;
imagesc(max_proj)
colormap(gray)
hold on
plot(squeeze(num_all(:,:,1)/conversion),squeeze(num_all(:,:,2)/conversion),'.')

if no_roi == 1;
    [x,y] = getpts;

    in = inpolygon(num_all(:,1,1),num_all(:,1,2),x*conversion,y*conversion);
    plot(squeeze(num_all(~in,:,1)/conversion),squeeze(num_all(~in,:,2)/conversion),'.b')
    plot(squeeze(num_all(in,:,1)/conversion),squeeze(num_all(in,:,2)/conversion),'.r')
elseif no_roi == 0;
else
    [x1,y1] = getpts;
    in1 = inpolygon(num_all(:,1,1),num_all(:,1,2),x1*conversion,y1*conversion);
    
    [x2,y2] = getpts;
    in2 = inpolygon(num_all(:,1,1),num_all(:,1,2),x2*conversion,y2*conversion);
    
    plot(squeeze(num_all(:,:,1)/conversion),squeeze(num_all(:,:,2)/conversion),'.b')
    plot(squeeze(num_all(in1,:,1)/conversion),squeeze(num_all(in1,:,2)/conversion),'.r')
    plot(squeeze(num_all(in2,:,1)/conversion),squeeze(num_all(in2,:,2)/conversion),'.g')
end

strflag = questdlg('Is/Are the ROI/s correct?');
if strcmp(strflag,'Cancel') == 1;
   exit
end
end

close(1)

%calculate displacements and speed
track = zeros(no_tracks,no_frames,9);
track(:,:,9) = num_all(:,:,4);
for i = 1:no_tracks;
    for j = 1:no_frames-1;
        track(i,j+1,1) = abs(num_all(i,j,1) - num_all(i,j+1,1));
        track(i,j+1,2) = abs(num_all(i,j,2) - num_all(i,j+1,2));
        track(i,j+1,3) = abs(num_all(i,j,3) - num_all(i,j+1,3));
        track(i,j+1,4) = sqrt(track(i,j+1,1)^2 + track(i,j+1,2)^2 + track(i,j+1,3)^2);
        track(i,j+1,5) = track(i,j+1,4)*(1/str2num(answer{4}));
        track(i,j+1,6) = sqrt((num_all(i,j+1,1) - num_all(i,1,1))^2 + (num_all(i,j+1,2) - num_all(i,1,2))^2 + (num_all(i,j+1,3) - num_all(i,1,3))^2);
        track(i,j+1,7) = track(i,j+1,4) + track(i,j,7);
        track(i,j+1,8) = track(i,j+1,6)/track(i,j+1,7);
    end
end

%********Align times*****************
for i = 1:no_tracks;
    if sum(isnan(track(i,:,1))) > 0;
        temp = squeeze(track(i,:,:));
        track(i,:,:) = NaN;
        track(i,min(temp(:,9))+1:max(temp(:,9))+1,:) = temp(1:find(isnan(temp(:,1)),1)-1,:);
    end
end

%process data ready for R plotting
info = zeros(no_tracks,no_frames,2);
concat = cat(3,info,num_all);
%number all cells by ID
for i = 1:no_tracks
    concat(i,:,1) = i;
end

%label all tracks with 0, then label ROIs (1 or 1&2)
concat(:,:,2) = 0;

if no_roi == 1;
    concat(in,:,2) = 1;
elseif no_roi == 0;
else
    concat(in1,:,2) = 1;
    concat(in2,:,2) = 2;
end


%normalise all morphological parameters
concat2 = concat;
for i = 1:no_tracks;
    if sum(~isnan(concat(i,:,20))) <= 2;  
    else
    for k = [10 11 12 20 21 22 23];
        for j = 1:no_frames;
        concat2(i,j,k) = concat(i,j,k)/concat(i,find(~isnan(concat(i,:,k)),1),k);
        end
    end
    end
end


%convert time to hours post injury
concat(:,:,6) = time_mat;

%Now for the msd/ad data
array2 = cat(3,concat(:,1:no_frames-1,6),msd_all);
array3 = cat(3,concat(:,1:no_frames-2,6),da_all);

%convert time to delta t
array2(:,:,1) = deltat_mat; 
array3(:,:,1) = deltat2_mat;

tracks2 = zeros(no_tracks,no_frames,4);
tracks2(:,1:no_frames-1,1:2) = array2;
tracks2(:,1:no_frames-2,3:4) = array3;

full = cat(3,concat,concat2(:,:,[10:12,20:23]),track(:,:,1:8),tracks2);

test = reshape(full,no_tracks*no_frames,42);

col_header1={'cellid','roi','x pos','y pos','z pos','time','x size','y size','z size',...
    'sph','round','conv','maxferet','ellipseA','ellipseB','ellipseC','yaw','pitch','roll','elongation','flatness',...
    'sa','vol','sphn','roundn','convn','elongn','flatnessn','san','voln'...
    'xdisp','ydisp','zdisp','inst speed','instspeed','disp','distance','directionality',...
    'msddeltat','msd','dadeltat','da'};     %Row cell array (for column labels)
fname1=[Folder_Name '_ggplot.csv']; %csv filename
writetable(cell2table([col_header1' num2cell(test')]'),fname1,'writevariablenames',0); %write data

test_mean = squeeze(mean(full,2,'omitnan'));

for i = [31:33,36:38];
   for j = 1:no_tracks;
    if sum(isnan(full(j,:,i))) >= no_frames - 2;%helps remove for any short tracks or tracks with little information 
        test_mean(j,i) = NaN; %CHECK THIS IS CORRECT
    else
        %now to find last non-nan value in each column which isn't all NaN
        A = full(j,:,i);
        B = ~isnan(full(j,:,i));
        test_mean(j,i) = A(find(B, 1, 'last'));
    end
   end
end

        
%finds the last non-NaN number in each 'last' parameter required e.g.x disp
%y disp z disp disp dist and directionality
for i = [31:33,36:38];
    A = full(:,:,i)';
    B = ~isnan(A);
 % indices
 Indices = arrayfun(@(x) find(B(:, x), 1, 'last'), 1:size(A, 2));
 % values
 Values = arrayfun(@(x,y) A(x,y), Indices, 1:size(A, 2));
 test_mean(:,i) = Values;
end

%calculating the DA AUC
for j = 1:no_tracks;
test_mean(j,42) = (1+full(j,1,42))*0.5*str2num(answer{4}); % initialise the AUC
for i = 1:no_frames-3;
    if full(j,i,42) ~= full(j,i,42);
        full(j,i,42) = 0;
    end
    if full(j,i+1,42) ~= full(j,i+1,42); %if NaN
        full(j,i+1,42) = 0;
    end
    test_mean(j,42) = test_mean(j,42) + 0.5*str2num(answer{4})*(full(j,i,42)+full(j,i+1,42));
end
end

col_header2={'roi','x size','y size','z size',...
    'sph','round','conv','maxferet','ellipseA','ellipseB','ellipseC','yaw','pitch','roll','elongation','flatness','sa','vol'...
    'sphn','roundn','convn','elong','flatnessn','san','voln'...
    'xdisp','ydisp','zdisp','inst speed','instspeed','disp','distance','directionality',...
    'msd','da AUC'};  
fname2 = [Folder_Name '_ggplot_mean.csv'];
writetable(cell2table([col_header2' num2cell(test_mean(:,[2,7:38,40,42])')]'),fname2,'writevariablenames',0)
   