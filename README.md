## PuntoMorph Readme - PLEASE READ!
High content analysis of phagocytic activity and cell morphology

## Installing the standalone executable on Windows
Download and run the PuntoMorph1.2.exe file to install PuntoMorph on your machine. This will require you to install the appropriate MATLAB runtime environment. Note that you will need a 64-bit Windows (and administrator access) to install the standalone software.<br />
Depending on the speed of your computer, PuntoMorph may take a few minutes to start the first time you run it. Please be patient. 

## Compiling standalone software for Mac or Linux
To compile a standalone executable for Mac and Linux, download the source code (including the MATLAB figure file) and all the m.file dependecies and compile on your local Mac or Linux machine. This can also be done for Windows. 

## Testing out your PuntoMorph installation
1- Download and unzip the sample image folder provided in the Github repository (demo_images.zip) <br />
2- Run PuntoMorph <br />
3- Point PuntoMorph's "Select Folder" button to select the folder containing the images (e.g. the folder you unzipped in #1) <br />
4- Select the appropriate settings for your image format and color channels <br />
5- Run PuntoMorph <br />
6- Examine the output Analysis.txt file using a file editor (e.g. Wordpad) or spreadsheet editor (e.g. Excel) <br />

## General usage considerations
1- Ensure that images are of reasonable quality and resolution <br />
2- Preferably include 30 or more valid cells per image to enable accurate exclusion of artifacts <br />
3- For phagocytic analysis, avoid excessive clumping of bead objects <br />
4- Examine the tracing overlay of analyzed images to identify potential errors <br />
5- Use the parameters reported in the results file to apply additional filters <br />
6- Every time you run the program on a folder, it will delete previous analysis and log files. If you want to save the output files from a previous run, move those files to a new folder before rerunning the program<br />

## Running PuntoMorph in MATLAB
If you plan to run PuntoMorph in MATLAB, download the PuntoMorph.m and PuntoMorph.fig files and place them in the path of your MATLAB installation. You will also need the following m-files (available from MATLAB File Exchange): <br />

circle_hough.m <br />
circle_houghpeaks.m <br />
circlepoints.m <br />
dlmcell.m <br />
imoverlap.m <br />
getfilenames.m  <br />

**p.s. We advise replacing all instances of "strcmp" in getfilenames.m with "strcmpi"**
