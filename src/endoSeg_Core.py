## Add some details about the code liscensing my name etc


################################################## IMPORT LIBRARIES ####################################################

# Style reference - functions have underscores, variables are camelCase


################################################### PLOTTING SETUP #####################################################
# This section of code is not super essential but it helps make figures more
# customizable.

# Define figure axis font sizes
plt.style.use('seaborn-muted')
SMALL_SIZE = 6
MEDIUM_SIZE = 8
BIGGER_SIZE = 10
plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


################################################ MAIN EndoSeg CLASSES ##################################################

class endoSegSegmentor:
    """This class handles the different types of segmentation that need to be
       performed by endoSeg. Mostly it is similar but depends on how the input
       masks are generated. Masks from ML will need a certain amount of cleanup.
       Masks that are hand drawn are likely good as is. dataType defines if the
       input data is multi channel or time series"""

    def __init__(self, outPath, path2Data, path2Mask, modelSource, dataType, dapiChannel):
        self.outDir = outPath
        self.modelSource = modelSource
        self.maskData = io.imread(path2Mask)
        self.data = np.squeeze(convert_to_Npy(path2Data))
        zz, yy, xx = self.data.shape
        self.zz = zz
        self.yy = yy
        self.xx = xx

        if dataType == 'multi channel':
            self.nucleusImage = self.data[dapiChannel, :, :]
            self.nucSeg = self.segment_nucleus()
        self.segResult = self.segment_from_mask()


    def segment_nucleus(self):
        """This function is used to segment DAPI labeled Nuclei from multichannel
           fluoresence data using simple morphological operations. The output is
           both the binary nucleus image (for removing nuclear signal) and the
           labelled nuclei."""

        thresh = threshold_otsu(self.nucleusImage)
        bw = closing(self.nucleusImage > thresh, square(4))
        cleared = remove_small_objects(clear_border(bw), 150)

        return cleared


    def segment_from_mask(self):
        """This function is used to clean up the predicted borders from machine learning
           or hand drawn borders and use them in segmentation."""

        maskData = rgb2gray(self.maskData)
        maskData = resize(maskData, (self.yy, self.xx))
        thresh = threshold_otsu(maskData)
        binMask = maskData > thresh

        if self.modelSource == 'handDrawn':
            borders = np.array(binMask)
            trimmedMaskData = remove_small_holes(binMask, 200)
            trimmedMaskData = remove_small_objects(trimmedMaskData, min_size=500)

        elif self.modelSource == 'ML':
            binMask = np.array(invert(1 * binMask))
            res = white_tophat(binMask, disk(5))
            binMask = dilation(binMask, disk(1))
            binMask = binMask - res
            borders = binMask
            trimmedMaskData = remove_small_holes(binMask, 200)
            trimmedMaskData = remove_small_objects(trimmedMaskData, min_size=500)

        cellLabels = ndi.label(trimmedMaskData)[0]
        wsSeg = watershed(trimmedMaskData, cellLabels, mask=borders)

        return wsSeg


    def create_overlay_image(self):
        """This function writes out an overlay image displaying the results of the
           segmentation with each cell labelled accordingly."""

        ## Get the location of all the cells and their ID and create the overlay image
        allCentx = []
        allCenty = []
        allIds =  []
        cellRegionProps = regionprops(self.segResult)
        for region in cellRegionProps:
            centY, centX = region.centroid
            allCentx.append(centX)
            allCenty.append(centY)
            allIds.append(str(self.segResult[int(centY, int(centX))]))

        imageLabelOverlay = label2rgb(self.segResult, self.data[3, :, :]*40, bg_label=0)

        ## Plot the overlay image and label the cells then save the image

        f, ax = plt.suplots(1, 1, figsize=(5, 5), dpi=300)
        ax.imshow(imageLabelOverlay, interpolation='none')
        for i, txt in enumerate(allIds):
            ax.annotate(txt, (allCentx[i], allCenty[i]), color='white', size=7)

        ax.set_title('segmentation')
        plt.tight_layout()
        plt.savefig(self.outDir + 'Segmentation_result.png')



# Rethink how to do this a little bit. Have the class now correspond to a single cell
# such that it gets passed the data from the segmentation class on a per cell basis.
# Can then do any analysis that needs doing on a cell by cell call from within a loop

## To do from 2022/02/10 - Need to write a function that actually integrates the signal within cells
##                       - Figure out what format to pass data to class to intially.
##                       - Which functions need to sit outside the class?


class endoSegAnalyser:
    """This class performs the different analysis steps of endoSeg on a per cell basis.
       This class can handle the two main data types from these experiments; either a
       single time point multiple channel image looking at localisation of different proteins,
       or a single channel, multiple time point dataset for analyzing calcium localisation
       and spikes. Input should be a masked cell."""

    def __init__(self, dataIn, mask, cellLabel, correctionType, segResult, channelsOfInterest,
                 dataType, outdir):
        self.data = dataIn  # Keep this as the data for convenience
        self.correctionType = correctionType
        self.segResult = segResult
        self.dataType = dataType  # Define whether time series of multichannel
        self.channelsOfInterest = channelsOfInterest  # This is important for multichannel image
        self.mask = mask
        self.cellLabel = cellLabel
        self.cellSignal = self.integrate_fluorescent_signal()
        self.correctedCellSignal = self.correct_decay()

    def calcOrientationMap(self):

        """This function takes in the  the segmentation and uses it to calculate the angle the cell is rotated
           away from the normal. Then, using the  centroid of the cell, it calculates a a distance map from the
           midpoint rotated according to the cells orientation which can be used to quantify protein concentrations
           at a given distance. cellToAnalyze is an integer."""

        # Retrieve the orientations (in radians) of the nuclei from region props

        cellRegionProps = regionprops(self.cellLabel)
        for region in cellRegionProps:
            centY, centX = region.centroid
            rotAngle = region.orientation

        print(centY, centX)

        if np.pi / 4 < rotAngle <= np.pi * (3 / 4):
            rotAngle -= np.pi / 2
        elif rotAngle > np.pi * (3 / 4):
            rotAngle -= np.pi

        rotAngle *= (180 / np.pi)
        print(str(rotAngle))

        # Calculate distance map rotated so that it is parallel to cell orientation
        yy, xx = np.indices((ySize, xSize))
        xx = xx - np.round(centX)
        yy = yy - np.round(centY)
        # This assumes flow is along the horizontal direction
        orientMap = xx  # rotate(yy, rotAngle, center=(centX, centY), order=0)

        return orientMap

    def getSectorInds(self, orientMap, sectorSize, downstreamEnd, upstreamEnd):

        """This funtion is used to define the indices of a sector or quadrant of the cell
           to perform analysis on. The input is the distance map calculated for a particular
           cell and indRange defines the regions of the map that need to be extracted"""
        downstreamEnd = int(downstreamEnd)
        upstreamEnd = int(upstreamEnd)

        sectorMap = np.zeros_like(orientMap)
        currentSector = 1
        for i in range(downstreamEnd, upstreamEnd + 1):
            if i != 0:
                if i <= (downstreamEnd + currentSector * sectorSize):
                    sectorMap[orientMap == i] = currentSector
                else:
                    currentSector += 1
                    sectorMap[orientMap == i] = currentSector

        return sectorMap

    def integrate_fluorescent_signal(self):
        """This will be the function that goes from an individual cell to a 1D fluoresence trace.
           Make sure this works for both data types, will probably need some case switching code."""


    def correct_decay(self):
        """Function for fitting a mixed exponential model to calcium signalling data to account for the
           signal loss due to fluoresence decay over time. If fit type is 1 the fit uses an asymmetric
           loss function to prevent over-fitting the curve to the data. This is to ensure that only the
           background is subtracted and spikes in calcium signal are preserved. If type is 2 the fit will
           use a fairly broad guassian kernal to estimate a smoothed background from the data directly."""
        if self.correctionType == 1:
            yHat = np.array(self.cellSignal)
            x = np.arange(0, len(self.cellSignal))

            params = Parameters()
            params.add('amp', value=20)
            params.add('decay', value=-1)
            params.add('a', value=0.01)
            params.add('b', value=-0.5)
            # params.add('m', value=-0.5)
            params.add('c', value=450)
            epsData = 1.0

            # do fit, here with leastsq model
            minner = Minimizer(residual, params, fcn_args=(x, yHat, epsData), nan_policy='omit')
            try:
                result = minner.minimize(method='leastsq')
                amp = result.params['amp']
                decay = result.params['decay']
                a = result.params['a']
                b = result.params['b']
                # m = result.params['m']
                c = result.params['c']
                # print(amp, decay, a, b, c)
                model = (amp * np.exp(x * decay)) + a * x ** 2 + b * x + c

            except ValueError:
                # Sometimes least squares fails to fit due to NaN. In those cases switch
                # to the Nelder-Mead
                result = minner.minimize(method='nelder')
                amp = result.params['amp']
                decay = result.params['decay']
                a = result.params['a']
                b = result.params['b']
                # m = result.params['m']
                c = result.params['c']
                # print(amp, decay, a, b, c)
                model = (amp * np.exp(x * decay)) + a * x ** 2 + b * x + c
                # model = x + c
                # model = np.median(cellSignal)
                print('fail')
            cellSignal = self.cellSignal - model
            cellSignal[cellSignal < 0] = 0

        elif self.correctionType == 2:
            # If least squares is consistently failing use a heavily smoothed curve to
            # estimate the background.
            model = smoothByConvolve(self.cellSignal, 31, 'hanning')
            cellSignal = self.cellSignal - model
            cellSignal[cellSignal < 0] = 0

        elif self.correctionType == 3:
            # Final option for background subtraction is to fit a polynomial to subsets
            # of the curve to calculate a contiguous background
            model = savitzky_golay(self.cellSignal, 51, 5)
            cellSignal = self.cellSignal - model
            cellSignal[cellSignal < 0] = 0

        return cellSignal




