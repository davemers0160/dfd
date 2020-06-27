
%% get the file

file_name = 'D:\IUPUI\Test_Data\Lidar\test_2018-08-09-19-30-22.bag';

% read in the file

fileID = fopen(file_name,'r','ieee-le');

% read in the string for the ROSBAG version
tline = fgetl(fileID);


%% read in a record

% 1. read in header length
header_length = fread(fileID,1,'uint32','ieee-le');
% 2. read in the header data
% 2.1 read in field length
field_length = fread(fileID,1,'uint32','ieee-le');

% 2.2 read in the field name value pair
field_nv = fread(fileID,field_length,'char*1','ieee-le');


% A = fread(fileID,sizeA,precision,skip,machinefmt)