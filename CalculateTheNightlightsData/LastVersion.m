
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reema AlHassan %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Feel free to contact me via my Email : Reema.amh97@gmail.com for any suggestion or problems about the code..
%%%% Supervisor : Hector G. LOPEZ-RUIZ 
%%%% Email: hector.lopez@kapsarc.org
%%%% This program count the "Sum of the lights" in the stalliate images %%%%%%%%%%%%%%
%%%% some functions are from the File exchange section in MathWorks you can find the Links below %%%%%
%%% before using it first make sure about the conditions below:
%%%  1. you have C compiler such as xcode.
%%%  2. Symbolic Math and Mapping toolboxes.
%%%  3. you put all the files of the zip folder in MATLAB folder
%%%  4. Run the file "insidepoly_install.m" in the zip folder
%%%  5. if you choosed choice 2(If you have the images already in a folder) make sure you didn't rename the Tiff images..
%%%  6. make sure about the shape folder it should contains ".shp" , ".shx" and ".dbf" files
%%%  7. the geometry of the .shp files should be "polygons"
%%% the unzipping may takes a long time because of the image size do not
%%% worry about that, however to speed it up make sure you have a good internet connection
%%% i tried to cover most of the cases which cuase errors ,please excuse any uncoverd..

%%%%%%%%%%%%%%%%%%%%%%%%%%% Reema AlHassan %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% inputs block start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Please choose :');%% ask the user to choose either to unzip the URLs or read the images directly from a folder
input1 = input('1: If you have a text file of URLs \n2: If you have the images already in a folder\n3: Exit \n');

   while(input1~=1)&&(input1~=2)&&(input1~=3) %% vildate user input the ..input should be 1 or 2 or 3   
         fprintf('\n');
 disp('Invalid choice your choice should be from 1-3');
disp('please choose:');
input1 = input('1: If you have a text file of URLs \n2: If you have the images already in a folder\n3: Exit \n');
   
   end  
   
   
  
%------------------------------- choice 1 block start --------------------------------%         
 if(input1==1)%% if the user choosed to unzip the URLs
     
       % it should be text file
disp('Select a text file please'); % taking the textfile name from the user 
[file,path] = uigetfile('*.txt', 'Select a text file please'); %to get the file from diffrent folders
fprintf('The text file is: ');
fprintf(file);
fprintf('\n');
while isequal(file,0)%% validate user input%%
%%if the user pressed cancel the file will be 0 and we need it to be real file
     disp('There are no text file..');
    disp('Select a text file please');
    [file,path] = uigetfile('*.txt', 'Select a text file please');   
end

textFullfile = fullfile(path,file); %% making the path 

%%% asking the user about the name of the folder where that the extracted file will be 
%%% it can be already existing or if it not the program will creat it  
folder = input('Enter the folder name (where the extracted files will be): ','s'); 
   end 
  
%------------------------------- choice 1 block end --------------------------------%   




  
%------------------------------- choice 2 block  start --------------------------------%   
         if(input1==2) %% if the user has the images already in a folder 
disp('Please select the folder which has the images'); 
TIFFFOLDER= uigetdir(); %% to get the folder path
TiffPattern = fullfile(TIFFFOLDER, '*.tif'); %% collect only the images with tif extension
tifFiles   = dir(TiffPattern);
         end 
%------------------------------- choice 2 block end --------------------------------%   


  

%------------------------------- choice 3 block start --------------------------------%   
if (input1==3) %% Exit the program by stop excution
   error('Program has been stopped, this error message was made on purpose to stop the excution')
end 
%------------------------------- choice 3 block end--------------------------------%   






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% in both choices 1 and 2 the user needs to specify the shape folder and the algorithm below %%%%%%

disp('Select the shape folder please'); 
shapeFolder= uigetdir(); %Shape folder
  
 
%%% which algorithm to apply in the intersection between the image and the shape
%%% "inpolygon" more precise in the most cases but if the shape contains large polygons it will take along time maybe more than 5 hours for each image  
%%% so if your polygons are large i suggest you to use "insidepoly" it is much faster in this case but it
%%% doesn't give the exact values as QGIS but close values ,however if your polygons are
%%% small i suggest you to use "inpolygon"...
disp('Please choose the algorithm number to get the "Sum of the lights":') 
fprintf('\n');
algorithm = input('1-inpolygon: more precise but it takes time for large polygons\n2-insidepoly:less precise but much faster in most of the cases\n');%%% validate user input


    while(algorithm~=1)&&(algorithm~=2) %% vaildate user input 
         fprintf('\n');
 
disp('Invalid input please enter either 1 or 2');
disp('please choose the algorithm number to get the "Sum of the lights":') 
algorithm = input('1:inpolygon\n2:insidepoly\n');
    end 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% inputs block end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% processing block start%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

%------------------------------- choice 1 block  start --------------------------------% 
   
if (input1==1) %% to complete unzipping for choice 1 
        
disp('Please wait the program is runnig ...');
str= fileread(textFullfile); % read the txt file line by line 
 
pattern = regexpi(str, ...
    ['((http|https|ftp|file)://|www\.|ftp\.)',...
    '[-A-Z0-9+&@#/%=~_|$?!:,.]*[A-Z0-9+&@#/%=~_|$]'], 'match'); % a pattern for reading urls
 

%%%%% unzipping the files in the output folder 

for k = 1:numel(pattern)
url=pattern{k};
ind=strfind(url,'.tif'); 
if(isempty(ind)==0) %% there is a .tif
lastword= url(ind(end)+1:end);
lastword=lower(lastword); %% making them both in a lower case 
%%%%%%%%%%%%%%%%%%%%% if the image is daily %%%%%%%%%%%%%%%%%%%%%
%%% this case only for the daily URLs since they don't contain ".tgz "
%%% to better understandig this is a sample of daily URL:
%%https://data.ngdc.noaa.gov/instruments/remote-sensing/passive/spectrometers-radiometers/imaging/viirs/mosaics//20180704/SVDNB_npp_d20180704.d.75N060W.rade9.tif 
if strcmp('tif',lastword)
    gunzip(url,folder);
end
 
else %%%%%%%%%%%%%%%%%%%%% else (yearly or monthly)%%%%%%%%%%%%%%%%%%%%%
%%% this case for the monthly and yearly because they are  ".tgz "
gunzip(url,folder);%% this one must to have 
untar(url,folder); %% this function unzip the .tgz
%%% to better understandig this is a sample of yearly URL:
%%%% https://data.ngdc.noaa.gov/instruments/remote-sensing/passive/spectrometers-radiometers/imaging/viirs/dnb_composites/v10//2015/SVDNB_npp_20150101-20151231_75N060W_v10_c201701311200.tgz

end

 end


%%%% extracting the TIFF images only for the daily, monthly and yearly
dirList =dir(folder);
names = {dirList.name}; % it contains all the files the TIFF images and others
  TIFFimages={}; % this will contain only the TIFF images
j=1;
for i=1:numel(names)
 
  if contains(names(i),'vcmcfg') && contains(names(i),'rade9h') && contains(names(i),'vcm-orm')==0  %%%% 1) checking for the monthly images
   
                TIFFimages(j)=names(i);
         j=j+1;
  end
  
  if  contains(names(i),'vcm-orm') && contains(names(i),'rade9')&& contains(names(i),'vcmcfg')==0  %%%% 2)  checking for the yearly images
 
              TIFFimages(j)=names(i);
        j= j+1;    
  end
  
   if contains(names(i),'vcm-orm')==0 && contains(names(i),'vcmcfg')==0 && contains(names(i),'rade9')==1 && contains(names(i),'vcm-ntl_v10')==0  && contains(names(i),'vcm_v10')==0 %%%% 3) checking for daily images

            TIFFimages(j)=names(i);
        j= j+1;    
   end
   %%added by me
   if contains(names(i),'web.stable')==0 && contains(names(i),'lights.avg')==0 %%%% 3) checking for dmsp images

            TIFFimages(j)=names(i);
        j= j+1;    
   end %%added by me
  
end  
 

%%%%%%%%%%%%% making a folder which contains only the TIFF images the 
%%%%%%%%%%%%% the folder name will be 'TiffImages+random number between 1-1000' 

done =0;
while(done==0)
rand = randi(1000); %% compute a random integer
rand=num2str(rand); %% convert it to string
TIFFFOLDER= strcat('TiffImges',rand); %% concatenate it with the word "TiffImges"
 if exist(TIFFFOLDER)==0 %% checking if the current folder has a folder with this name (TiffImages+random number) 
names1TiffImagesFolder=mkdir(TIFFFOLDER); %% if not then the folder will be created 
done=1; %% to exit from the loop
break
 else %% if the folder has this name then  this random number will be skipped 
     continue
 end
end

%%% moving the Tiff images to the created folder
for z=1:numel(TIFFimages)
 fullPathOfTiff=fullfile(folder,TIFFimages{z});
m= movefile(fullPathOfTiff,TIFFFOLDER); 
end
 
fprintf('TIFF images have been collected in the folder: ');%%% displying a message for the user 
fprintf(TIFFFOLDER);
fprintf('\n');

TiffPattern = fullfile(TIFFFOLDER, '*.tif');
tifFiles   = dir(TiffPattern);
    
    end      
%------------------------------- choice 1 block end --------------------------------% 
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% the code below for both choices 1 and 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% here to make sure that the shape folder contains the (shp , shx, dpf ) files

shp=fullfile(shapeFolder,'*.shp');%% finding .shp files
shx=fullfile(shapeFolder,'*.shx');%% finding .shx files
dbf=fullfile(shapeFolder,'*.dbf');%% finding .dbf files

shpFiles=dir(shp);
shxFiles=dir(shx);
dbfFiles=dir(dbf);

%%% if one of these (shp , shx, dpf ) files not inside the folder
   while((length(shpFiles)==0)||(length(shxFiles)==0)||(length(dbfFiles)==0)) 
        
        fprintf('\n');
        fprintf('\n');
        
disp('The folder does''t contain one of these extensions : .shp ,.shx and .dbf');
disp('Please make sure it has these files and select it again');
shapeFolder= uigetdir(); %Shape folder
      
       
shp=fullfile(shapeFolder,'*.shp');%% finding .shp files
shx=fullfile(shapeFolder,'*.shx');%% finding .shx files
dbf=fullfile(shapeFolder,'*.dbf');%% finding .dbf files

shpFiles=dir(shp);
shxFiles=dir(shx);
dbfFiles=dir(dbf);

   end
   
%%% This loop go through all the .shp files to find the one which has a polygon Geometry
%%% this loop doesn't take all the .shp files in a folder it only
%%% checks where is the correct .shp one (which has the polygons only not points or lines)
%%% so if you have multiple .shp files put them please in separate folders.
 for j = 1:length(shpFiles)
     
baseShape=shpFiles(j).name;
shape = shaperead(fullfile(shapeFolder,baseShape));
geo=shape(1).Geometry; %% here i checked only for the first one becuse it's impossible that the file has mex geometry so checking for the first one is enough
geo=lower(geo);

if strcmp(geo,'polygon')==1  %% here to make sure that the geometry is polygon and nothing else
break

else    
continue
end

 end
   
%%%%%%%%%% this loop will go through all the images and map them shape to all of them %%%%%%%%%%%%%%
for k = 1:length(tifFiles)
  baseFileName = tifFiles(k).name;
  tifFullFileName = fullfile(TIFFFOLDER, baseFileName);

   
numOfpolygons=size(shape,1); %% the number of rows represent the number of polygons 
total=0; 

a = zeros(numOfpolygons,1); %% 'a' is a matrix which will contain the sum of the lights for each polygon
%%% this loop will go through all the polygons and in the shape file and map them with the image
for i=1:numOfpolygons 

try
    
%%%% clipping the image to a square each time based on the x and y coordinates of the shape 
%%%%% this step to prevent the program from crashing by minimizing the image size 
[values,imageX,imageY,info] =geoimread(tifFullFileName,shape(i).X,shape(i).Y);

shapeX = shape(i).X(1:end-1);%% taking the x coordinates and exclude the NaN from the matrix 
shapeY = shape(i).Y(1:end-1);%%taking the y coordinates and exclude the NaN from the matrix 
 
[gridsX,gridsY] = meshgrid(imageX,imageY); %% convert the image coordinates to grids befor intersect it with the shape coordinates 



if(algorithm==1) %% if the user choosed algorthem 1 
mask  =inpolygon(gridsX,gridsY,shapeX,shapeY);  %% this function intersect the image grids with the shape coordinates 
%%% then it makes a "Mask" of the intersection area the mask represent binary matrix of zeros and ones
end 


if(algorithm==2)%% if the user choosed algorthem 2 
mask  =insidepoly(gridsX,gridsY,shapeX,shapeY); %% this function intersect the image grids with the shape coordinates 
%%% then it makes a "Mask" of the intersection area the mask represent binary matrix of zeros and ones
end 


result= mask.*values; %% multiply the mask with the image values  element by element 
sumOfValues= sum(result(:));%% sum all the values of the resulting mstrix to get the sum of the lights
sumOfValues=double(sumOfValues)
a(i) = sumOfValues;
total=total+sumOfValues; %% after each iteration sum the value with the previous values to get the total of the whole shape
total=double(total);


catch
%%% if the erea of the polygon too small for the function to specify it 
sumOfValues= NaN %% then here it will be replaced with a "NaN" 
a(i) = sumOfValues;
total=total+0;
total=double(total);

   continue 
end

end
 
disp('The total:');
total

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% processing block end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output block start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %%% here to cut the date from the images to make the syntax name for the
 %%% csv file ....
 %%% for example here are some names of tiff images
 %%% the yearly: SVDNB_npp_20150101-20151231_75N060W_vcm-orm_v10_c201701311200.avg_rade9.
 %%% the monthly:SVDNB_npp_20170701-20170731_75N060W_vcmcfg_v10_c201708061230.avg_rade9h.tif
 %%% the daily SVDNB_npp_d20180704.d.75N060W.rade9.tif
 %%% so the file name will be based on the kind of the image
 
 %------------------------- Yearly images output syntax ---------------------------------------%

 if contains(baseFileName,'vcm-orm')==1 %% becuase each yearly image has this word

if contains(baseFileName,'vcm-orm-ntl')==1 %%% because there are 2 types of the yearly images 
  type='vcm-orm-ntl';                      %%% and here we are identify them in the output
  type1='ntl';
else 
    type='vcm-orm';
    type1='';
    
end
    
SumOfLightColumn = extractBetween(baseFileName,'SVDNB_npp_','-','Boundaries','exclusive');

SumOfLightColumn=strcat('Sum',SumOfLightColumn);
SumOfLightColumn=char(SumOfLightColumn);
SumOfLightColumn = SumOfLightColumn(1:end-4); %% cutting only the year number so the output here will be for example Sum2015

SumOfLightColumn1=SumOfLightColumn; %% temp

SumOfLightColumn=cellstr(SumOfLightColumn);
SumOfLightColumn=strcat(SumOfLightColumn,type1);
SumOfLightColumn=char(SumOfLightColumn);



 end
 
 
 %------------------------- Monthly images output syntax ---------------------------------------%
 
 if contains(baseFileName,'vcmcfg')==1
  
SumOfLightColumn = extractBetween(baseFileName,'SVDNB_npp_','-','Boundaries','exclusive');
SumOfLightColumn=strcat('Sum',SumOfLightColumn);
SumOfLightColumn=char(SumOfLightColumn);
SumOfLightColumn = SumOfLightColumn(1:end-2);%% the output here will be for example Sum201702 where 2017 is the year and 02 is the month
  
type='vcmcfg';

SumOfLightColumn1=SumOfLightColumn; %% temp
  
 end
  
 
 %------------------------- Daily images output syntax ---------------------------------------%
 
    %%% because here each daily image has ".d." only
 if contains(baseFileName,'.d.')==1 && contains(baseFileName,'vcmcfg')==0 && contains(baseFileName,'vcmcfg')==0
    
SumOfLightColumn = extractBetween(baseFileName,'SVDNB_npp_d','.','Boundaries','exclusive');
SumOfLightColumn=strcat('Sum',SumOfLightColumn);
SumOfLightColumn=char(SumOfLightColumn);%% the output here will be for example Sum2018704 2
       type='';  

SumOfLightColumn1=SumOfLightColumn; %% temp      
      end
 
shape = shaperead(fullfile(shapeFolder,baseShape)); %% here i read the shape again to get the shape info after the previous iteration
shapeToTable=struct2table(shape,'AsArray',true); %% here setting AsArray to true to allow diffrent number of rows
n=num2cell(a); %% convert the a from double array to cell array 
s=setfield(shapeToTable,SumOfLightColumn,n); %% "n" array is "a" array but after convertion to a cell array 
w = table2struct(s); %% convert the table to struct 
shapefile1=fullfile(shapeFolder,baseShape);%% taking the full path of the shape 
shapewrite(w,shapefile1); %% writing the info to the shape again

s.X = []; %%% removing the x coordinates from the table in the csv file only not from the shape
s.Y = []; %%% removing the y coordinates from the table in the csv file only not from the shape
s.BoundingBox = []; %%% removing the BoundingBox from the table in the csv file only not from the shape


csvfFileName=cellstr(SumOfLightColumn1);
csvfFileName=strcat(SumOfLightColumn1,type);

csvfFileName=strcat(csvfFileName,baseShape,'.csv'); %% making the syntax for the file name 
csvfFileName=char(csvfFileName);%% convert it to character because it is a cell array
writetable(s,csvfFileName);%% write the table in csv file


end 
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output block end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%% here all the functions which i have re-used from MathWorks %%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function 1 geoimread %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%i use it because it is  more efficient with large Geotiff images(satellite images)
%%%% the link : https://ww2.mathworks.cn/matlabcentral/fileexchange/46904-geoimread....
%%% Author Info: 
% 
% (c) Aslak Grinsted 2014 (http://www.glaciology.net/)
%     & Chad A. Greene (http://chadagreene.com/)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function 1 geoimread %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,x,y,I]=geoimread(filename,varargin)
 
%% Set defaults: 
 
usegeocoordinates = false; 
returnNONsubsetImage = true; 
buffer_x = 0; 
buffer_y = 0; 
 
%% Input checks: 
 
% Check for mapping toolbox: 
% assert(license('test','map_toolbox')==1,'geoimread requires Matlab''s Mapping Toolbox.')
 
% Check file type: 
assert(isnumeric(filename)==0,'Input filename must be a string.') 
[~,~,ext] = fileparts(filename);
switch upper(ext)
    case {'.JP2' '.JPEG2000' '.GEOJP2'}
        I = jp2tiffinfo(filename);
    case {'.TIF' '.TIFF' '.GTIF' '.GTIFF'}
       
       I = robustgeotiffinfo(filename);
    otherwise
        error('Unrecognized image file type. Must be tif, tiff, gtif, gtiff, jp2, jpeg2000, or geojp2.')
end
 
 
% Parse optional inputs: 
if nargin>1
    returnNONsubsetImage = false; 
    assert(nargin>2,'If you specify an xlim or latlim, you must specify a corresponding ylim or lonlim.')
    
    % Parse limits: 
    xlimOrLatlim = varargin{1}(:);
    ylimOrLonlim = varargin{2}(:);
    
    assert(isnumeric(xlimOrLatlim)==1,'Input xlim or latlim must be numeric.')
    assert(isnumeric(ylimOrLonlim)==1,'Input ylim or lonlim must be numeric.')
    
    % Assume geo coordinates if no input limits exceed normal lat/lon values: 
    if max(abs(xlimOrLatlim))<=90 && max(abs(ylimOrLonlim))<=360
        usegeocoordinates = true; 
    end
    
    % Parse buffer: 
    if nargin>3
        buffer_m = varargin{3}; 
        assert(isnumeric(buffer_m)==1,'Buffer value must be either a scalar or a two-element array.')
        assert(numel(buffer_m)<3,'Buffer value must be either a scalar or a two-element array.')
        buffer_x = buffer_m(1); 
        if isscalar(buffer_m)
            buffer_y = buffer_m(1);
        else
            buffer_y = buffer_m(2); 
        end
    end
    
    if nargin>4
        error('Too many inputs in geoimread.') 
    end    
    
end
 
 
%% Begin work: 
 
% Get pixel coordinates of full (non-subset) image: 
[x,y]=robustpixcenters(I);
 
% Set xlim and ylim depending on user inputs: 
if returnNONsubsetImage
    xlimOrLatlim = x(:);
    ylimOrLonlim = y(:); 
end
 
if usegeocoordinates
    % lat/lon limits switch to x/y limits here: 
    if ~strcmp(I.ModelType,'ModelTypeGeographic')
        assert(license('test','map_toolbox')==1,'Mapping toolbox needed to project between lat/lon limits and x,y limits. Specify limits in x,y coordinates.')
        [xlimOrLatlim,ylimOrLonlim]=projfwd(I,xlimOrLatlim,ylimOrLonlim);
    end
end

xlim = [min(xlimOrLatlim)-buffer_x max(xlimOrLatlim)+buffer_x];
ylim = [min(ylimOrLonlim)-buffer_y max(ylimOrLonlim)+buffer_y];
 
 
% Rows and columns of pixels to read: 

rows=find((y>=ylim(1))&(y<=ylim(2)));
cols=find((x>=xlim(1))&(x<=xlim(2)));


%% Display messages if region of interest is partly or wholly outside the image region:
 
if xlim(1)<min(x)||xlim(2)>max(x)
    disp('geoimread limits extend beyond the available image output in the x direction.')
end
 
if ylim(1)<min(y)||ylim(2)>max(y)
    disp('geoimread limits extend beyond the available image output in the y direction.')
end
 
if isempty(rows)||isempty(cols)
    error('No image coordinates can be found inside the specified limits.')
    
end
 
%% Load region of interest: 
reductionlevel=0;
if reductionlevel==0
    rows=sort(rows([1 end]));
    cols=sort(cols([1 end]));
else
    %% Load region of interest:
    dpx=2^reductionlevel;
    rows=round(rows/dpx);cols=round(cols/dpx);
    
    rows=sort(rows([1 end]));
    cols=sort(cols([1 end]));
end
x=x(cols(1):cols(end));
y=y(rows(1):rows(end));
 
A=imread(filename,'PixelRegion',{rows cols});
 
 
%% Update info structure to more accurately reflect the new image: 
 
if nargout == 4
    I.FileSize = numel(A); 
    I.Height = size(A,1); 
    I.Width = size(A,2); 
    try
        I.TiePoints.WorldPoints.X = x(1);
        I.TiePoints.WorldPoints.Y = y(1);
        I.SpatialRef.RasterSize = [size(A,1),size(A,2)];
        I.RefMatrix(3,1) = x(1);
        I.RefMatrix(3,2) = y(1);
        I.BoundingBox = [min(x) min(y); max(x) max(y)];
        I.CornerCoords.X = [min(x) max(x) max(x) min(x)];
        I.CornerCoords.Y = [max(y) max(y) min(y) min(y)];
        %TODO: check whether GTRasterTypeGeoKey is RasterPixelIsArea or RasterPixelIsPoint
        I.CornerCoords.Row = .5 + [0 0 size(A,1) size(A,1)]; %TODO: is this .5 always true?  
        I.CornerCoords.Col = .5 + [0 size(A,2) size(A,2) 0];
        [I.CornerCoords.Lat,I.CornerCoords.Lon] = projinv(I,I.CornerCoords.X,I.CornerCoords.Y);
        I.GeoTIFFTags.ModelTiepointTag(4) = x(1);
        I.GeoTIFFTags.ModelTiepointTag(5) = y(1);
        I.SpatialRef.XLimWorld = [min(x),max(x)];
        I.SpatialRef.YLimWorld = [min(y),max(y)];
    catch,end
end
 

%% Clean up: 
 
if nargout==0
    imshow(A,'XData',x,'YData',y)
    axis xy
    clear A x y I
end
 
 


 
function I=jp2tiffinfo(fname)
% code to extract embedded geotiff from jp2 file...
 
fid=fopen(fname,'r','ieee-be');
%JPg2000 info: http://www.jpeg.org/public/15444-1annexi.pdf
%geojp2 info: http://www.lizardtech.com/download/geo/geotiff_box.txt
 
while ~feof(fid)
    lbox=fread(fid,1,'uint32');
    type=fread(fid,[1 4],'uint8=>char');
    lbox=lbox-8;
    if lbox==1
        lbox=fread(fid,1,'uint64');lbox=lbox-8;
    end
    if strcmp(type,'uuid')
        uuid=fread(fid,[1 16],'uint8'); lbox=lbox-16;
        geo=[177 75 248 189 8 61 75 67 165 174 140 215 213 166 206 3];sprintf('\xb1\x4b\xf8\xbd\x08\x3d\x4b\x43\xa5\xae\x8c\xd7\xd5\xa6\xce\x03');
        if all(uuid==geo)
            foutname=fullfile(tempdir,'tifembeddedinjp2.tif');
            fout=fopen(foutname,'w');
            contents=fread(fid,lbox,'uint8');lbox=0;
            fwrite(fout,contents);
            fclose(fout);
            fclose(fid);
            I=robustgeotiffinfo(foutname);
            m=imfinfo(fname); % a little silly to use imfinfo when i already have a tag reader
            I.Height=m.Height;
            I.Width=m.Width;
            delete(foutname);
            return
        end
    end
    fseek(fid,lbox,0);
end
fclose(fid);
 
%--- BELOW is to make it more robust if no mapping toolbox ---
 


end 

function I = robustgeotiffinfo(fname)
if license('test','map_toolbox')
    I=geotiffinfo(fname);
else
    I=imfinfo(fname);

end
end
function [x,y]=robustpixcenters(I)
if license('test','map_toolbox')
    [x,y]=pixcenters(I);
else
    %I have not read documentation... but this only works for rectilinear systems.
    assert(I.ModelPixelScaleTag(3)==0,'unexpected ModelPixelScaleTag format.');
    assert(all(I.ModelTiepointTag(1:3)==0),'unexpected ModelTiepointTag format.');
    x=((0:I.Width-1)-I.ModelTiepointTag(1))*I.ModelPixelScaleTag(1)+I.ModelTiepointTag(4);
    y=((0:I.Height-1)-I.ModelTiepointTag(2))*-I.ModelPixelScaleTag(2)+I.ModelTiepointTag(5);
end

end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function 2 insidepoly%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% it is for the algorithm 2 
%%% it was bulit using MEX files.
%%% before using it first make sure about the conditions below:
%%%  1. you have C compiler. 
%%%  2. Mathmatics toolbox and Mapping toolbox.
%%%  3. you put all the files of the zip folder(in the link) in MATLAB folder
%%% except the file which has the function because i have put it here already
%%% Function link : https://ww2.mathworks.cn/matlabcentral/fileexchange/27840-2d-polygon-interior-detection?s_tid=answers_rc2-1_p4_MLT
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%     Original: 06-Jun-2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function 2 insidepoly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [inpoly onboundary] = insidepoly(varargin)

% Default options
options = struct('tol', 'auto', ...
                 'presortflag', 'auto');
                 
optionloc = cellfun('isclass',varargin,'char');
if any(optionloc)
    optionloc = find(optionloc,1,'first');
    mainargin = varargin(1:optionloc-1);
    % retreive property/value pairs
    for k=optionloc:2:nargin-1
        options.(varargin{k})= varargin{k+1};
    end
else
    mainargin = varargin;
end

edges = NaN;
if length(mainargin)<=3
    [xy P] = deal(varargin{1:2});
    x = xy(:,1);
    y = xy(:,2);
    if length(mainargin)==3 % edge provided
        edges = varargin{3};
    end
    Px = P(:,1);
    Py = P(:,2);
else
    [x y Px Py] = deal(varargin{1:4});
     if length(mainargin)==5 % edge provided
        edges = varargin{5};
     end
end

edgeprovided = ~isnan(edges);

% Orginal size of the data
sz = size(x);

% Reshape in columns
Px = Px(:);
Py = Py(:);
x = x(:);
y = y(:);

% Cast to double of one of them is (required by Mex)
isxdbl = isa(x,'double');
isydbl = isa(y,'double');
isPxdbl = isa(Px,'double');
isPydbl = isa(Py,'double');
doubleengine = isxdbl || isydbl || isPxdbl|| isPydbl;
if doubleengine
    if ~isxdbl, x = double(x); end
    if ~isydbl, y = double(y); end
    if ~isPxdbl, Px = double(Px); end
    if ~isPydbl, Py = double(Py); end   
end

Pxmin = min(Px);
Pxmax = max(Px);
Pymin = min(Py);
Pymax = max(Py);

% Select the tolerance for determining on-boundary points
if isequal(options.tol,'auto')
    ontol = 1e-9*max(Pxmax-Pxmin,Pymax-Pymin);
else
    ontol = options.tol;
end

if ~edgeprovided
    % We don't want duplicate end vertices
    if Px(1)==Px(end) && Py(1)==Py(end)
        Px(end) = [];
        Py(end) = [];
    end
    % Linear wrap around the points
    Px1 = Px;
    Py1 = Py;
    Px2 = Px([2:end 1]);
    Py2 = Py([2:end 1]);
else
    Px1 = Px(edges(:,1),:);
    Py1 = Py(edges(:,1),:);
    Px2 = Px(edges(:,2),:);
    Py2 = Py(edges(:,2),:);
end

% Filter data outside the rectangular box
inpoly = x>=Pxmin-ontol & x<=Pxmax+ontol & ...
         y>=Pymin-ontol & y<=Pymax+ontol;
x = x(inpoly);
y = y(inpoly);

if isequal(options.presortflag,'auto')
    n = size(x,1);
    m = size(Px1,1);
    % We don't do presorting for small size polygon or small size data
    % Empirical law from experimental tests (Bruno)
    presortflag = (m>32) && (n*m > 25000);
else
    presortflag = options.presortflag;
end

if presortflag
    % Sort the array in x, and find the brackets. This is used to find
    % easily which data points have abscissa fall into the abcissa bracket
    % of each edge of the polygon.
    [x ix first last] = presort(Px1, Px2, x);
    y = y(ix); % arrange y in the same order
    
    % Call mex engine
    %
    if doubleengine
        [in on] = insidepoly_dblengine(x, y, Px1, Py1, Px2, Py2, ontol, ...
                                    first, last);
                               
    else % single arrays
       [in on] = insidepoly_sglengine(x, y, Px1, Py1, Px2, Py2, ontol, ...
                                        first, last);
                                  
    end
    
    % Restore the original order
    in(ix) = in;
    on(ix) = on;
else % No presorting
    % Call mex engine, without presorting
    if doubleengine
       [in on] = insidepoly_dblengine(x, y, Px1, Py1, Px2, Py2, ontol);
        
    else % single arrays
        [in on] = insidepoly_sglengine(x, y, Px1, Py1, Px2, Py2, ontol);
         
    end
end

in = in | on;

if nargout>=2
    onboundary = inpoly;
    onboundary(inpoly) = on;
    % Reshape to original size
    onboundary = reshape(onboundary, sz);
end

inpoly(inpoly) = in;
% Reshape to original size
inpoly = reshape(inpoly, sz);

end % insidepoly

%%
function [xsorted ix first last] = presort(Px1, Px2, x)
% (Px1, Px2) abscissa of vertices, x abscissa of data, they are supposed
% to be ranged in column.
% Return:
%   xsorted as sort(x) = x(ix)
%   "first" & "last" indexes such that, for each vertice point
%       min(Px1,Px2) <= xsorted(first:last) <= max(Px1,Px2)

% left and right brackets of the segment
Pmin = min(Px1,Px2);
Pmax = max(Px1,Px2);

nvertices = size(Px1,1);

% We seek to see how x interveaves with Pmin by sorting the ensemble
[trash is] = sort([Pmin; x],1); %#ok
isdata = is>nvertices; % tail index, i.e., belong to data abscissa x
anchor = find(~isdata);
% Get the sorted data alone
ix = is(isdata)-nvertices;
xsorted = x(ix); % sorted x
% Index of the first element in xsorted such that
%   xsorted(first)>=sort(Pmin)
first = anchor-(0:nvertices-1).';
% Rearrange first corresponds to the original order
ip = is(anchor);
first(ip) = first;

% determine how Pmax interleaves with xsorted, i.e.,
% index of the last element in xsorted such that xsorted(first)<=Pmax
% Note: in case of draw in binning edges, HISTC must return the last edge
[trash last] = histc(Pmax, [xsorted; inf]); %#ok

end % presort



    

  
