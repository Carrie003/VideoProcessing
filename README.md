# VideoProcessing 

## File Structure
    src/ - Directory with source code
      Makefile - Unix/Mac makefile for building the project with "make". 
      imagepro.[vcproj/sln/suo] - Project file for Visual Studio 2005 on Windows
      imgpro.cpp - Main program, parses the command line arguments, and calls the appropriate image functions
      R2Image.[cpp/h] - Image class with processing functions (this is the only file that you need to edit)
      R2Pixel.[cpp/h] - Pixel class 
      R2/ - A library of useful 2D geometric primitives
      jpeg/ - A library for reading/writing JPEG files
    input/ - Contains example input images. 
    output/ - Es empty to start -- it will contain the images produced by your program (see below)
    runme.bat - a script (for Windows) that you will fill in to demonstrate execution of your program
    runme.sh - same as <code>runme.bat, but for Mac OS X
    
## Compilation
### Windows
If you are developing on a Windows machine and have Visual Studio
installed, use the provided project solution file (assn1.sln) in the
src/ directory to build the program. 
### Mac/Linux machine
If you are developing on a Mac or Linux machine, cd into the src/ directory and type "make". In either case, an executable called imgpro (or imgpro.exe) will be created in the src/ directory.


## How To Run
```
# Clone this repository
$ git clone https://github.com/XinyuYang/VideoProcessing.git

# Go into the repository
$ cd VideoProcessing

# To see all the options of operations
$ src/imgpro -help

# Image processing (Drag the image into the input folder, put the method name in the end of the command line such as brighntess,blur,harris...)
$ src/imgpro input/testpattern.jpg output/testpattern_brighntess_0.5.jpg -brightness 0.5

# Video processing (Drag the images of the video into the videoinput folder; create a new folder called "videooutput" to save the output sequence)
$ src/imgpro -video

```

## Run your own video
### Prerequisites
FFmpeg
```
# The FFmpeg command line syntax to convert a video into a JPG sequence
$ ffmpeg -i video.mp4 -qscale:v 1 videoinput/input%07d.jpg

# The FFmpeg command line syntax to convert a JPG sequence into a video (24fps)
$ ffmpeg -r 24 -f image2 -s 1920x1080 -i videooutput/output%07d.jpg -vcodec libx264 -crf 25 -pix_fmt yuv420p test.mp4
```

## Demo
![Alt text](img-demo/skyreplacement.png?raw=true "Title")
Link to all test videos:      
https://drive.google.com/drive/folders/14GVWiPm0vweJ_soODPyko9QYaXaEwfuM?usp=sharing

## Authors
* Xinyu Yang
* Weiqiu You
* Megan Zhao
