# endoSeg
EndoSeg is a set of functions for calculating different image statistics from fluorescent images of endothelial cells. It is set up to handle either time course images (where fluoresence is being monitored in a single channel as the series progresses) or multi-channel images (where fluorescent is being monitored and compared in different fluorescent channels). EndoSeg can handle a variety of image formats although it was written with .czi and .tif files in mind. Dependecies can be found in the associated .yml files
Author: Marcus Gallagher-Jones  
Institution: The Rosalind Franklin Institute  
Email: marcusgj13@gmail.com

## Running endoSeg

The easiest way to interact with endoSeg is via the provided jupyter notebook. Analysis in endoSeg is split into three distinct sections: Data loading, preprocessing and segmentation and feature extraction and analysis. EndoSeg requires two inputs. The first is a stack of fluoresecent images (preferrably a 3D data cube where the first dimension is the size of the image stack and the second and third dimensions are the individual image dimensions (i.e. ZXY)). The second is an image that defines the boundaries of the cells in the field of view. It can be either hand drawn or derived from another segmentation software. The outline image does not need to be a binary image as endoSeg should be able to handle the conversion.



