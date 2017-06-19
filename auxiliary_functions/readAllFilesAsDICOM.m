function  status = readAllFilesAsDICOM(directory_name,subFolder)
% Read all DICOM-files in specified folder, sorts and writes all studies
% and series to .mat files.
% INPUT: 
%       - directory_name: path of folder to search
%       - subFolder: 0 scan only this folder (default), 1 scan subfolders
%       as well
% OUTPUT: Found studies and series are written in 'directory_name'
%       - STRUCT: data
%           data.data: Scan volume in double
%           data.info: DICOM header
%           data.SliceLocation: [x, y, z] in scanner coordiante system
%
% VERSION: 0.91
%
% Copyright Anders Nymark Christensen, anym@dtu.dk, DTU Compute 20140217

% Add neccesary utilities to path
addpath('C:\Dropbox\ArbejdeDTU\MATLAB_utilities')

% Check System
if ispc
    slash = '\';
    septoken = ';';
elseif isunix
    slash = '/';
    septoken = ':';
else
    error('System unknown. Try hardcoding slash and septoken')
end

% check first input
if ~strcmp(directory_name(end),slash)
    directory_name = [directory_name,slash];
end

% Check second input
if nargin > 1
    if subFolder ~= 0 && subFolder ~=1
        error('second input shoul be 0 or 1 if specified')
    end
else
    subFolder = 0;
end

%Change to chosen folder
cd(directory_name)
disp(['Current folder is: ',directory_name])

%% Find alle studies and series

% Get filelist
FileList = getAllFiles(directory_name,subFolder);

% If no files abort
if isempty(FileList)
    status = 'No files in folder';
    return
end

% Allocate memory
Modality = cell(length(FileList),1);
SeriesID = cell(length(FileList),1);
StudyID = cell(length(FileList),1);
ImPosPatient = zeros(length(FileList),3);

% Find all studies and series
for i=1:size(FileList,1)
    filename = FileList{i,1};
    try
        fileinfo = dicominfo(filename);
        IM = dicomread(filename);
        if isempty(IM)
           error('')
        end
        fileinfo = dicominfo(FileList{i,1});
        Modality{i,1} = fileinfo.Modality;
        SeriesID{i,1} = fileinfo.SeriesInstanceUID;
        StudyID{i,1} = fileinfo.StudyInstanceUID;
        if isfield(fileinfo,'ImagePositionPatient')
            ImPosPatient(i,:) = fileinfo.ImagePositionPatient;
        else
            ImPosPatient(i,:) = [nan nan nan];
        end
    catch
        Modality{i,1} = '';
        SeriesID{i,1} = '';
        StudyID{i,1} = '';
        ImPosPatient(i,1) = nan;
    end
end

% If no DICOM files abort
if isempty(Modality)
    status = 'No DICOM files in folder';
    return
end

[Cstudy,IAstudy,ICstudy] = unique(StudyID);
[Cseries,IAseries,ICseries] = unique(SeriesID);
[Cmodality,IAmodality,ICmodality] = unique(Modality);

% IM = zeros(fileinfo.Columns,fileinfo.Rows,size(FileList,1));
% disp(['Study name is: ',fileinfo.PatientID])

status = 'Alas my Captain: No images here!';
%% Read and write all

% Loop over all studies and all series
for indexStudy = 1:length(Cstudy)
    if ~isempty(Cstudy{indexStudy})
        for indexSeries = 1:length(Cseries)
            if ~isempty(Cseries{indexSeries})
                
                % All files in series
                seriesFiles = FileList(ICseries == indexSeries);
                % Allocate memory for files
                fileinfo = dicominfo(seriesFiles{1,1});
                
                if strcmpi(fileinfo.Modality , 'PT') % PET
                    disp('PET')
                    if ~isfield(fileinfo,'NumberOfTimeSlices')
                        disp('Static...')
                    
                 
                        % Get data
                        for indexFiles =  1:length(seriesFiles)
                            if isfield(fileinfo,'RescaleSlope') && isfield(fileinfo,'RescaleIntercept')
                                IM(:,:,indexFiles) = fileinfo.RescaleSlope .*...
                                    double(dicomread(seriesFiles{indexFiles,1})) + fileinfo.RescaleIntercept;
                            else
                                IM(:,:,indexFiles) = double(dicomread(seriesFiles{indexFiles,1}));
                            end
                            fileinfo = dicominfo(seriesFiles{indexFiles,1});
                            
                            times(indexFiles,1)=str2double(fileinfo.StudyTime(5:end))+...
                                60*str2double(fileinfo.StudyTime(3:4))+...
                                3600*str2double(fileinfo.StudyTime(1:2));
                            
                            times(indexFiles,2)=str2double(fileinfo.SeriesTime(5:end))+...
                                60*str2double(fileinfo.SeriesTime(3:4))+...
                                3600*str2double(fileinfo.SeriesTime(1:2));
                            
                            times(indexFiles,3)=str2double(fileinfo.AcquisitionTime(5:end))+...
                                60*str2double(fileinfo.AcquisitionTime(3:4))+...
                                3600*str2double(fileinfo.AcquisitionTime(1:2));
                            
                            times(indexFiles,4)=str2double(fileinfo.ContentTime(5:end))+...
                                60*str2double(fileinfo.ContentTime(3:4))+...
                                3600*str2double(fileinfo.ContentTime(1:2));
                            
                            times(indexFiles,5) = fileinfo.ActualFrameDuration / 1000;
                            
                        end
                        
                        SliceLocation = ImPosPatient(ICseries == indexSeries,:);
                        [~, ind] = sort(SliceLocation(:,3),'ascend');
                        SliceLocation = SliceLocation(ind,:);
                        IM = IM(:,:,ind);
                        times = times(ind,:);
                        
                        data.data = IM;
                        data.info = fileinfo;
                        data.times = times;
                        data.SliceLocation = SliceLocation;
                        save(['Study',num2str(indexStudy-1),'_Series',num2str(indexSeries-1),'_',fileinfo.Modality,'_',fileinfo.SeriesDescription,'.mat'],'data','-v7.3')
                    else
                        disp('Dynamic...')
                        for indexFiles =  1:length(seriesFiles)
                            IM(:,:,indexFiles) = squeeze(fileinfo.RescaleSlope *...
                                double(dicomread(seriesFiles{indexFiles,1})) + fileinfo.RescaleIntercept);
                            
                            fileinfo = dicominfo(seriesFiles{indexFiles,1});
                            
                            times(indexFiles,1)=str2double(fileinfo.StudyTime(5:end))+...
                                60*str2double(fileinfo.StudyTime(3:4))+...
                                3600*str2double(fileinfo.StudyTime(1:2));
                            
                            times(indexFiles,2)=str2double(fileinfo.SeriesTime(5:end))+...
                                60*str2double(fileinfo.SeriesTime(3:4))+...
                                3600*str2double(fileinfo.SeriesTime(1:2));
                            
                            times(indexFiles,3)=str2double(fileinfo.AcquisitionTime(5:end))+...
                                60*str2double(fileinfo.AcquisitionTime(3:4))+...
                                3600*str2double(fileinfo.AcquisitionTime(1:2));
                            
                            times(indexFiles,4)=str2double(fileinfo.ContentTime(5:end))+...
                                60*str2double(fileinfo.ContentTime(3:4))+...
                                3600*str2double(fileinfo.ContentTime(1:2));
                            
                            times(indexFiles,5) = fileinfo.ActualFrameDuration / 1000;
                        end
                        
                        
                        SliceLocation = ImPosPatient(ICseries == indexSeries,:);

                        
%                             [~, ind] = sort(times(:,4),'ascend');
%                             times = times(ind,:);
%                             for i = 1:length(ind)
%                                 IM2{i,1} = IM{ind(i),1};
%                             end
                            [~, ind] = sort(times(:,4),'ascend');
                            times = times(ind,:);
                            
                            uniqueTimes = unique(times(:,5));
                            for i = 1:length(uniqueTimes)
                                SliceLocationIM{i,1} = SliceLocation(times(:,5) == uniqueTimes(i),:);
                                [~, ind] = sort(SliceLocationIM{i,1}(:,3),'ascend');
                                dummy = IM(:,:,times(:,5) == uniqueTimes(i));
                                IM2{i,1} = dummy(:,:,ind);
                            end
                            

                        [~, ind] = sort(SliceLocation(:,3),'ascend');
                        SliceLocation = SliceLocation(ind,:);
                        IM = IM(:,:,ind);
                            
                            data.data = IM2;
                            data.info = fileinfo;
                            data.times = times;
                            data.SliceLocation = SliceLocation;
                            clear IM2
                            
                            save(['Study',num2str(indexStudy-1),'_Series',...
                                num2str(indexSeries-1),'_',fileinfo.Modality,'_',...
                                fileinfo.SeriesDescription,'.mat'],'data','-v7.3')
                            status = 'My Captain: Images here!';
                    end
                    
                else % CT MR
                    IM = zeros(fileinfo.Height,fileinfo.Width,length(seriesFiles));
                    disp('CT / MRI')
                    disp(['Series size: ',num2str(fileinfo.Height),'  ',num2str(fileinfo.Width),'  ',num2str(length(seriesFiles))])
                    % Get data
                    
                    for indexFiles =  1:length(seriesFiles)
                        if isfield(fileinfo,'RescaleSlope') && isfield(fileinfo,'RescaleIntercept')
                            IM(:,:,indexFiles) = fileinfo.RescaleSlope *...
                                double(dicomread(seriesFiles{indexFiles,1})) + fileinfo.RescaleIntercept;
                        else
                            IM(:,:,indexFiles) = double(dicomread(seriesFiles{indexFiles,1}));
                        end
                    end
                    SliceLocation = ImPosPatient(ICseries == indexSeries,:);
                    [~, ind] = sort(SliceLocation(:,3),'ascend');
                    SliceLocation = SliceLocation(ind,:);
                    IM = IM(:,:,ind);
                    
                    
                    data.data = IM;
                    data.info = dicominfo(seriesFiles{indexFiles,1});
                    data.SliceLocation = SliceLocation;
                    save(['Study',num2str(indexStudy-1),'_Series',num2str(indexSeries-1),'_',fileinfo.Modality,'_',fileinfo.SeriesDescription,'.mat'],'data','-v7.3')
                    status = 'My Captain: Images here!';
                end
            end
            clear IM fileinfo ind data times
        end
    end
end



%% TEST

% n=1;
% for i=1:length(FileList)
%     filename = FileList{i,1};
%      if strcmpi(filename(end-2:end),'dcm') || strcmpi(filename(end-2:end),'ima') % Check if DICOM-file
%         fileinfo = dicominfo(FileList{i,1});
%         if isfield(fileinfo,'ActualFrameDuration') %isfield(fileinfo,'NumberOfFrames')
%      
%             test(n,1) = fileinfo.ActualFrameDuration;
%             n = n+1;
%         end
%      end
% end
% 
% for i = 1:length(test)
%     data = dicomread(FileList{test(i),1});
%     size(data)
% end