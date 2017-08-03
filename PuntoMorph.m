function varargout = PuntoMorph(varargin)
% PuntoMorph MATLAB code for PuntoMorph.fig
%      PuntoMorph, by itself, creates a new PuntoMorph or raises the existing
%      singleton*.
%
%      H = PuntoMorph returns the handle to a new PuntoMorph or the handle to
%      the existing singleton*.
%
%      PuntoMorph('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PuntoMorph.M with the given input arguments.
%
%      PuntoMorph('Property','Value',...) creates a new PuntoMorph or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PuntoMorph_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PuntoMorph_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PuntoMorph

% Last Modified by GUIDE v2.5 09-Nov-2016 15:06:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PuntoMorph_OpeningFcn, ...
    'gui_OutputFcn',  @PuntoMorph_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before PuntoMorph is made visible.
function PuntoMorph_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PuntoMorph (see VARARGIN)

% Choose default command line output for PuntoMorph
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PuntoMorph wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = PuntoMorph_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    %open logfile
    fid = fopen('logFile','a+');
    
    set(handles.pushbutton1,'String','Running...','Enable','off');
    set(handles.popupmenu4,'Enable','off');
    set(handles.popupmenu6,'Enable','off');
    set(handles.popupmenu7,'Enable','off');
    set(handles.popupmenu8,'Enable','off');
    set(handles.selectfolder,'Enable','off');
    set(handles.edit2,'Enable','off');
    
    pic=get(handles.popupmenu4,'value');
    chan1=get(handles.popupmenu6,'value');
    chan2=get(handles.popupmenu7,'value');
    chan3=get(handles.popupmenu8,'value');
    folder=get(handles.edit2,'String');
    
    if(pic==2)
        picformat='.jpg';
    elseif(pic==3)
        picformat='.tif';
    elseif(pic==4)
        picformat='.png';
    elseif(pic==5)
        picformat='.bmp';
    end
    
    channels=[chan1,chan2,chan3];
    
    if (pic==1 || min(channels)==1 || length(unique(channels))< 3 || strcmp(folder,'Please provide path to folder containing images'))
        set(handles.edit6,'String','Please make appropriate selections and try again');
    else
        %%
        warning('off', 'Images:initSize:adjustingMag');
        
        perc=1; %Optional: set this to ignore a given percentage of the upper and lower normally distributed data from analysis
        cell_nuc_ratio=1.25; %minimum ratio of accepted cell size to mean nucleus size. Gets rid of cell-like artifacts
        
        cd(folder);
        
        %Intialize results matrix
        result(1,1:12)={'Folder' 'Image' 'Cell number' 'Part of normal distribution' 'Area' 'Convex Area' 'Splay Index' 'Perimeter' 'Convex Perimiter' 'Eccentricity' 'Ruffling' 'Beads per cell'};
        
        %Find all images and their paths under current folder
        dnfiles = getfilenames(pwd,horzcat('*',picformat));  %Specify directory here
        [Missingmdb,~] = strtok(dnfiles,'.'); %Extracts first part of filename
        DBFile(1,1:length(dnfiles)) = {''}; %Initialize vars for speed
        DBPath(1,1:length(dnfiles)) = {''};
        DBFolder(1,1:length(dnfiles))={''};
        
        for z = 1:1:length(dnfiles)
            string = char(Missingmdb(1,z));
            string = fliplr(string);
            [DBFile(1,z), DBPath(1,z)] = strtok({string},'\');
            [DBFolder(1,z), DBPath(1,z)] = strtok(DBPath(1,z),'\');
            DBFile(1,z) = {fliplr(char(DBFile(1,z)))};
            DBPath(1,z) = {fliplr(char(DBPath(1,z)))};
            DBFolder(1,z) = {fliplr(char(DBFolder(1,z)))};
        end
        
        for filecounter = 1:1:length(dnfiles)
            
            %navigate to folder being currently analyzed
            cd (horzcat(DBPath{1,filecounter},DBFolder{1,filecounter}));
            
            %create folder for traced images
            Analysis=exist('Trace');
            if(Analysis == 7)
            else
                mkdir Trace;
            end
            
            set(handles.edit6,'String',(horzcat('Analysing ',DBFile{1,filecounter},' in ',DBFolder{1,filecounter})));
            fprintf(fid,'%s\n',horzcat('*** Analysing ',DBFile{1,filecounter},' in ',DBFolder{1,filecounter},'...'));
            fprintf(fid,'%s\n','Preparing image for analysis...');
            
            inputfile=horzcat(DBFile{1,filecounter},picformat);
            
            %read image and separate into its RGB channels
            I=imread(inputfile);
            
            axes(handles.axes2)
            imshow(I);
            pause(0.1);
            
            %Set the blue channel to either nuclei (B), cellbodies
            %(R), or beads (G)
            if(chan1==2)
                B = I(:,:,3);
            elseif(chan1==3)
                R = I(:,:,3);
            elseif(chan1==4)
                G = I(:,:,3);
            end
            
            %Set the red channel to either nuclei (B), cellbodies
            %(R), or beads (G)
            if(chan2==2)
                B = I(:,:,1);
            elseif(chan2==3)
                R = I(:,:,1);
            elseif(chan2==4)
                G = I(:,:,1);
            end
            
            %Set the green channel to either nuclei (B), cellbodies
            %(R), or beads (G)
            if(chan3==2)
                B = I(:,:,2);
            elseif(chan3==3)
                R = I(:,:,2);
            elseif(chan3==4)
                G = I(:,:,2);
            end
            
            %enhance contrast and threshold image
            R(R<(max(max(R))/10))=0; %delete any pixels that have less than 5% of the max intensity
            R_eq = adapthisteq(R,'cliplimit',0.03,'Distribution','rayleigh','Alpha',0.7);
            bw = im2bw(R_eq, graythresh(R_eq));
            
            %fill holes inside image objects
            bw1 = imfill(bw,'holes');
            %fudge image to connect small breaks
            bw2 = bwdist(bw1) <= 2;
            %smooth contours
            bw3 = imopen(bw2, ones(2,2));
            
            %threshold nucleus channel and pick nuclei
            bw5 = im2bw(B, graythresh(B));
            nuclei = imfill(bw5,'holes');
            nuclei = bwdist(nuclei) <= 2; %fudge nuclei a little
            nuclei = imopen(nuclei, ones(2,2)); %smooth nuclei contours
            nuclei = bwareaopen(nuclei, 100); %clean out tiny objects/specs
            
            fprintf(fid,'%s\n','Identifying nuclear centers...');
            
            %fit sizes to 2-normal distribution and guess mean area for single
            %nucleus
            
            %get rid of nuclei touching border
            nuclei=imclearborder(nuclei);
            
            nucs=bwlabel(nuclei);
            nuc_labels=unique(nucs);
            nuc_labels(1)=[]; %get rid of zero
            nuc_sizes(1:numel(nuc_labels))=0; %initialize nuclear size vector
            
            %if nuclei are detected, execute the rest of the code.
            %Otherwise skip this image and move to the next one.
            
            if(isempty(nuc_labels))
                
                fprintf(fid,'%s\n',horzcat('No cells detected. Skipping image...'));
                
            else
                
                for k = 1:numel(nuc_labels)
                    nuc_sizes(k)=length(nucs(nucs==k));
                end
                
                nuc_sizedist=gmdistribution.fit(nuc_sizes',2,'Regularize', 1e-5);
                
                %mean area of single nucleus
                nuc_avgsize=min(nuc_sizedist.mu);
                
                %%find individual nuclei
                mean_nucradius=round(sqrt(nuc_avgsize/pi));
                radii=(round(mean_nucradius*0.5)+1):1:round(mean_nucradius*1.5); %create range of radii roughly between half mean radius (ensuring it's > 0) and 1.5 mean radius
                h = circle_hough(nuclei, radii, 'same', 'normalise');
                peaks = circle_houghpeaks(h, radii, 'nhoodxy', (4*mean_nucradius+1), 'nhoodr', (4*mean_nucradius+1));
                
                nuc_centers=nuclei;
                nuc_centers(:)=0;
                nuc_locs=nuc_centers;
                
                %draw a cross over nuclear centers and overlay on image
                for peak = peaks
                    nuc_centers(peak(2), peak(1))=1;
                    nuc_locs(peak(2)-2:1:peak(2)+2,peak(1))=1;
                    nuc_locs(peak(2),peak(1)-2:1:peak(1)+2)=1;
                end
                
                %erode nuclei and trace them
                se = strel('disk',2,0);
                nuc_zone=imerode(nuclei,se);
                nuc_perim=bwperim(nuc_zone); %mark nuclear perimeters
                
                %apply watershed algorithm to segment cells by the nuclei within them
                R_eq_c = imcomplement(R);
                I_mod = imimposemin(R_eq_c, ~bw3 | nuc_centers);
                cells = watershed(I_mod);
                %ignore cells touching border
                cells=imclearborder(cells);
                %fill holes in cell watershed segments
                cells=imfill(cells,'holes');
                %remove tiny objects smaller than 1.5x mean nucleus size
                cells=bwareaopen(cells,round(cell_nuc_ratio*(nuc_avgsize)));
                
                %get rid of cells that don't have a nucleus inside them
                Label= bwlabel(cells); %count cells
                
                for x=1:1:(length(unique(Label))-1)
                    nuc_inside=sum(nuclei(Label==x));
                    if(nuc_inside < (round(nuc_avgsize))/2) %delete cells with less than 1/2 mean nucleus content
                        cells(Label==x)=0;
                    end
                end
                
                %draw cell perimiters
                cells_perim=bwperim(cells);
                
                %Recount cells and get their regionprops
                Label = bwlabel(cells);
                s = regionprops(cells, 'Orientation', 'MajorAxisLength', 'ConvexArea', 'Area','MinorAxisLength', 'Eccentricity', 'ConvexImage', 'Centroid', 'Perimeter', 'Image');
                
                %get intensities and sizes of all cells
                
                %initialize parameters
                intensities(1:numel(s))=0;
                sizes(1:numel(s))=0;
                cellnums(1:numel(s))=0;
                
                for k = 1:numel(s)
                    %get average intensity in cell body channel
                    intensities(k)=median(R(Label==k));
                    sizes(k)=length(Label(Label==k));
                    cellnums(k)=k;
                end
                
                fprintf(fid,'%s\n','Identifying valid cells...');
                
                %%fit intensities to a normal distribution and find the lowest and
                %%highest 1 percentile boundaries
                intensitiesdist=fitdist(intensities','Normal');
                xlowperc= intensitiesdist.mu-tinv(0.995,length(intensities))*(intensitiesdist.sigma);
                xhighperc= intensitiesdist.mu+tinv(0.995,length(intensities))*(intensitiesdist.sigma);
                
                %%fit intensities to a normal distribution and find the lowest and
                %%highest 1 percentile boundaries
                sizesdist=fitdist(sizes','Normal');
                xlowperc2= sizesdist.mu-tinv(0.995,length(sizes))*(sizesdist.sigma);
                xhighperc2= sizesdist.mu+tinv(0.995,length(sizes))*(sizesdist.sigma);
                
                %%eliminate outliers
                
                exclusions=unique(horzcat(cellnums(intensities<xlowperc),cellnums(intensities>xhighperc),cellnums(sizes<xlowperc2),cellnums(sizes>xhighperc2)));
                
                %find beads within cells
                Gpix=im2bw(G); %convert G to bw
                
                if(max(max(Gpix)>0))
                    fprintf(fid,'%s\n','Detected beads, analysing...');
                    
                    G_thres=G;
                    G_thres(G<(max(G(:))/5))=0;  %get rid of pixels with less < 20% max intensity (I have found this to be important because of green background outside cells)
                    bw_beads = im2bw(G);
                    
                    s_beads = regionprops(bw_beads,'Area','Centroid');
                    bead_labels=bwlabel(bw_beads);
                    bintensities(1:length(s_beads))=0;
                    
                    %find intensities of bead objects in the image
                    for b=1:1:length(s_beads)
                        
                        text(s_beads(b).Centroid(1,1)+3,s_beads(b).Centroid(1,2)+3,num2str(b),'FontSize',14,'Color','red');
                        %get the sum of bead intensities in the area corresponding to this
                        %object in the original image
                        bintensities(b)=sum(G_thres(bead_labels==b));
                        
                    end
                    
                    %fit intensities to 2-normal distribution and guess mean intensity for
                    %single bead
                    try
                        beadist=gmdistribution.fit(bintensities',2);
                        %mean area of single bead
                        bavgint=min(beadist.mu);
                    catch
                        bavgint=0; %if we can't find the average intensity of a single bead, avoid analyzing this image and generating wrong data
                        fprintf(fid,'%s\n',horzcat('Could not determine bead properties. Skipping image...'));
                    end
                    
                    %estimate number of beads per bead object
                    beads_present(1,length(s_beads))=0;
                    
                    for b=1:1:length(s_beads)
                        %divide total intensity for each object by mean intensity of single
                        %bead
                        beads_present(b)=round(bintensities(b)/bavgint);
                    end
                    
                end
                
                if(bavgint~=0)%do not analyze if we could not determine bead properties
                    
                    fprintf(fid,'%s\n','Annotating image...');
                    %overlay perimiters, hulls, and bead centers over cells
                    
                    overlay = imoverlay(I, cells_perim, [0 1 0]);
                    overlay = imoverlay(overlay, nuc_perim, [0 0 1]);
                    overlay = imoverlay(overlay, nuc_locs);
                    hulls=bwconvhull(cells,'Objects',8);
                    hulls_perim=bwperim(hulls);
                    overlay = imoverlay(overlay, hulls_perim);
                    
                    %write cell numbers on image
                    
                    for k = 1:numel(s)
                        c = s(k).Centroid;
                        overlay=insertText(overlay,[c(1),c(2)],k,'TextColor','white','BoxOpacity',0,'FontSize',18);
                    end
                    
                    %save figure
                    imwrite(overlay,horzcat(DBPath{1,filecounter},DBFolder{1,filecounter},'\Trace\','Trace_',strrep(inputfile,'.','-'),'.png'));
                    
                    %begin analyzing cell parameters
                    
                    fprintf(fid,'%s\n','Writing analysis output...');
                    
                    for k = 1:numel(s)
                        
                        ConvexImage=s(k).ConvexImage;
                        convex_perim=bwperim(ConvexImage);
                        convex_perimeter=sum(sum(convex_perim));
                        
                        %count bead centers from Hough transform method (beadnum) and also from intensity method (beadsincell)
                        if(max(max(Gpix)>0))
                            
                            %using intensity method
                            beadsincell=unique(bead_labels(Label==k));
                            beadsincell=beadsincell(beadsincell>0);
                            
                            totalbeads=0; %initialize parameter
                            
                            for h=1:1:length(beadsincell)
                                totalbeads=totalbeads+beads_present(beadsincell(h)); %sum up total number of guessed beads in each beady object
                            end
                        else
                            totalbeads=0;
                        end
                        
                        %markup whether script recommends inclusion or exclusion of cell
                        if(ismember(k,exclusions))
                            part_of=0;
                        else
                            part_of=1;
                        end
                        
                        %calculate Ruffling
                        clearvars perim_path perim_smoothed;
                        
                        smoothing_window=2*mean_nucradius; %let's just use the diameter of the nucleus for this
                        
                        ag=s(k).Image;
                        
                        ag=bwperim(ag);
                        
                        [x, y]=find(ag);
                        boundary=bwtraceboundary(ag,[x(1),y(1)],'N');
                        
                        perim_path(1:length(boundary))=0;
                        
                        sp=regionprops(s(k).Image,'Centroid');
                        centroid=sp.Centroid;
                        
                        for z=1:1:length(boundary)
                            perim_pixel=[boundary(z,1),boundary(z,2)];
                            X=[centroid;perim_pixel];
                            perim_path(z)=pdist(X,'euclidean');
                        end
                        
                        perim_path=perim_path-perim_path(1);
                        
                        Cterm=perim_path(end-(smoothing_window-2):end);
                        perim_path_flanked=horzcat(Cterm,perim_path);
                        
                        perim_smoothed=tsmovavg(perim_path_flanked','s',smoothing_window,1);
                        perim_smoothed(1:smoothing_window-1)=[];
                        perim_smoothed=circshift(perim_smoothed,1-smoothing_window/2);
                        
                        Ruffling=((s(k).Perimeter-convex_perimeter)/convex_perimeter)*(s(k).MinorAxisLength/s(k).MajorAxisLength);
                        
                        result=vertcat(result,horzcat(DBFolder(1,filecounter),DBFile(1,filecounter),(num2cell([k, part_of, s(k).Area, s(k).ConvexArea, (s(k).ConvexArea-s(k).Area)/s(k).ConvexArea, s(k).Perimeter, convex_perimeter, s(k).Eccentricity, Ruffling, totalbeads]))));
                    end
                end
                
                clearvars -except R G B fid handles dnfiles DBFile DBPath DBFolder perc result cell_nuc_ratio folder picformat chan1 chan2 chan3
                
            end
        end
        
        cd(folder);
        fprintf(fid,'%s\n','Done!');
        
        if(isempty(dnfiles))
            set(handles.edit6,'String','No images found');
        else
            dlmcell(horzcat('Analysis',regexprep(num2str(fix(clock)),'[^\w'']',''),'.txt'),result);
            set(handles.edit6,'String','Analysis Complete');
            set(handles.edit2,'String','Please provide path to folder containing images');
        end
    end
    
catch err
    
    set(handles.edit6,'String','Error: please check all settings and try again');
    set(handles.edit6,'String',DBFile{1,filecounter});
    set(handles.edit2,'String','Please provide path to folder containing images');
    
    cd(folder);
    %open file
    fprintf(fid,'%s\n',err.message);
    
    % following lines: stack
    for e=1:length(err.stack)
        fprintf(fid,'%sin %s at %i \n',err.stack(e).name,err.stack(e).line);
    end
    
    % close logfile
    fclose('all');
    
end

blackimage(1:90,1:100)=0;
axes(handles.axes2)
imshow(blackimage);

% close logfile
fclose('all');

set(handles.pushbutton1,'String','Run Analysis');
set(handles.pushbutton1,'Enable','on');
set(handles.popupmenu4,'Enable','on');
set(handles.popupmenu6,'Enable','on');
set(handles.popupmenu7,'Enable','on');
set(handles.popupmenu8,'Enable','on');
set(handles.selectfolder,'Enable','on');
set(handles.edit2,'Enable','on');

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in selectfolder.
function selectfolder_Callback(hObject, eventdata, handles)
% hObject    handle to selectfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder=uigetdir('C:\');
set(handles.edit2,'String',folder);

cd(folder);

%Clear old analysis file
oldanalysis=getfilenames('./','*Analysis*');
if(~isempty(oldanalysis))
    delete(oldanalysis{1});
end

%Clear old logfile
oldlogfile=getfilenames('./','*logfile*');
if(~isempty(oldlogfile))
    delete(oldlogfile{1});
end

guidata(hObject, handles);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over selectfolder.
function selectfolder_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to selectfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on key press with focus on selectfolder and none of its controls.
function selectfolder_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to selectfolder (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6

% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7

% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8

% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
