# matlab-mastrostack
Ma(e)stroStack: a Matlab class to automatically align and stack astro-photography images
 
 ![Image of MastroStack](https://github.com/farhi/matlab-mastrostack/blob/master/doc/mastrostack.png)
 
 Purpose
 -------
 
   This class gets a list of images, and automatically determines bright stars as control points. These are followed along pictures, and used to build an affine transformation at constant scale (e.g. a rotation and translation). All images are then stacked. The images can be given as file names (may include wildcards), or matrices from e.g. imread, and support both RGB and gray images. As stars are used for the alignment, this method is suited for deep sky images, but not for planetary imaging.
   
   This function does not make use of the phase correlation technique, but operates directly on the images. It assumes that at least two bright stars remain on each picture, but tolerate translations and rotations, as well as stars/control points that appear/disappear on the field of view. The procedure also ignores very sharp peaks, assuming they originate from dead pixels (and would then be static over images, leading to traces).
   
   It is highly recommended to specify a **dark** frame filename, which will be subtracted to all images, in order to e.g. remove most hot/dead pixels. To get such a 'dark', use your camera alone and shoot with the cap on, same settings as the other pictures (duration, ISO). Image should be full black.

  You may as well specify a **flat** frame filename, which will be divided to all images, in order to counter vignetting (darker borders). To get such a 'flat', shoot with the scope pointing at a uniform area (sky - blue or cloudy, white wall). Adapt the duration so that you get a rather gray image (not full black or white).
 
 Syntax
 -----------
  
**ma = mastrostack;**

  start the user interface
  
**ma = mastrostack(images)**

  loads images without setting their type
  
**ma = mastrostack(light, dark)**

**ma = mastrostack(light, dark, flat)**

  loads light, dark (background) and flat (scope response) images, and label them.
      
 Importing images
 ----------------
  
  Start with:
  
  ```matlab
    ma=mastrostack;
  ```
   
   Then press the **Return** key on the main interface. A **Drop Files Here** button appears in the lower left side. Drag and drop your Dark, Flat and Light images there. Images having 'dark' or 'flat' in their path/file name are marked as such automatically. You may alternatively use the File menu items.
   
   Supported image formats include JPG, PNG, TIFF, FITS. 
   
   If you have installed **readraw** (see Credits section below), you may as well directly import RAW camera images. This is highly recommended, as it retains much more information from the camera than the generated JPEG images, which proves to be essential for subtracting the Dark image (background), and revealing faint objects.
   
   Alternatively, you may use DCRAW on each RAW image with command
  
  ```
  dcraw -T -4 -t 0 <file.RAW>
  ```
   
   To do things quickly, you can just select the Compute/Stack menu item. This will set the Reference frame, compute the master Dark and Flat frames, align and stack images. Then go to the Stacking section below (you may skip the next section).
   
 Preparing the Stacking
 ----------------------
    
  The following steps are all optional.

  After importing the files, you should label them using the 'Image/Mark as...' menu items. You can navigate within images with the Image/Goto menu item, and the arrow keys, or the mouse wheel. 'Bad' images can be skipped (ignored). To use them back, set their type to 'light'. You should then compute the master Dark and Flat images (Compute menu).

  It is recommended to zoom onto specific features (e.g. a set of stars) to check visually for their sharpness. De-select the Zoom tool, and scan through images using the left arrow key, and press the 'I' key to mark images to be ignored, such as those blurred. To reset the plot, press the Return key.

  You can select the Reference image, which will be used as template for stacking. If not defined, the first image in the list will be used as such when stacking.

  Optionally, use the Compute/Align item to compute the images control points (stars) which also corrects for background and scope response when dark and flat are defined. This procedure computes a metric to quantify the sharpness. 
   
 Stacking
 --------
 
   When ready, use the Compute/Stack menu item. If the Alignment has not been executed previously, it is achieved for each image. The final image is then shown and written to disk. Use e.g. Lightroom, RawTherapee, DarkTable to enhance contrast.
   
 Improving the sharpness and contrast
 ------------------------------------

After the first Align or Stack procedure, the sharpness has been computed for all images. It is then possible to plot it in order to identify images of lower quality.

You may then select the 'Compute/Select on sharpness' menu item which metric is highest for clearer images, and low for blurred/moved ones. The axis is the image index.

You can use the following actions on that plot window. The mouse buttons define a rectangle selection.

- **LEFT**       button   selected images are set to SKIP/IGNORE
- **RIGHT**      button   selected images are set to LIGHT
- **SHIFT-LEFT** button   open first image in selection
- **LEFT/UP**    arrow    go to previous image
- **DOWN/RIGHT** arrow    go to next image
- **G**          key      goto image (selection from dialogue)
- **I/S**        key      current image is set to IGNORE/SKIP
- **L**          key      current image is set to LIGHT
- **X/Q/ESC**    key      abort 

You may as well check visually for the best images either from the main window (e.g. left/right arrows and press I key when the image is bad), or from the sharpness plot window with the same keys. It is recommended to zoom onto a portion of the image, to follow a few stars in shape. Then mark all blurred/strange star spots as Ignore/Skip.
   
 Notes
 -----
 
   If you have difficulties in stacking (some images do not have enough control points), relax e.g. the translation tolerance, using the menu item 'Compute/Set tolerances'. A typical value is 0.01 (1%), but higher values may help (e.g. 0.02). You can also increase the number of control points, and define the typical area of dead/hot pixels: control points/stars with lower extent area (usually 5-9) will be ignored for the definition of control points. A value of 0 will not handle dead pixels. These options usually require to restart the 'Align' procedure before Stacking.
   
   In case the main interface is closed, get it back with 
   
  ```matlab
    plot(ma)
  ```
  
  To save your session, with image references, labels, dark, flat, ignored images, and any stacked result, use:
  
  ```matlab
  save my_session
  ```
  
  which you can load back and display with:
  
  ```matlab
  load my_session
  plot(ma);
  ```
  
 
 Using commands (scripting)
 --------------------------
 
  Using commands allow to script/automate the procedure.
 
 ```matlab
    % create Ma(e)stroStack and import images
    ma=mastrostack('path/to/images/*.JPG','path/to/darks/*.JPG','path/to/flats/*.JPG');
 
    % stack. The first 'light' image will be used as Reference for stacking
    stack(ma);
```
 
 Methods
 -------
 
  - **about**           display a dialogue box
  - **correct**         correct and image for dark (background) and flat (vignetting)
  - **cpselect**        ALIGN images on reference
  - **diff**            compute difference of an image with reference
  - **delete**          delete the mastrostack, and clear memory
  - **exist**           check if an image has already been loaded
  - **getdark**         compute the master Dark
  - **getflat**         compute the master Flat
  - **help**            open the help page
  - **imread**          load an image and return its matrix and information
  - **label**           label an image as light, dark, flat or skip
  - **listdlg**         display a selector for images
  - **load**            load an image and return its information
  - **mastrostack**     create the Ma(e)strostack session
  - **plot**            plot the user interface and images
  - **stack**           STACK light images, correcting with dark and flat.
  
 Installation:
 -------------

  You only need Matlab, no external toolbox.
  
  Copy the directory and navigate to it. Then type from the Matlab prompt:

  ```matlab
  addpath(pwd)
  ```
  
  If you also have **readraw** installed and available (see below), you will be able to import
  RAW camera images. 
  
 Credits
 -------
  
  - **ReadRaw** a Matlab RAW camera reader <https://github.com/farhi/matlab-readraw>
  - **astrotnstack** automatically align and stack astro-photography pictures <https://github.com/farhi/astrotnstack>. This is the initial prototype of Ma(e)stroStack.
  - optimal rotation/translation between points <http://nghiaho.com/?page_id=671>
  - dndcontrol <https://fr.mathworks.com/matlabcentral/fileexchange/53511-drag---drop-functionality-for-java-gui-components>
  - blurMetric <https://fr.mathworks.com/matlabcentral/fileexchange/24676-image-blur-metric>

  See also:

  - LxNstack https://sites.google.com/site/lxnstack/home
  - Deep Sky Stacker http://deepskystacker.free.fr/french/
  - Rot'n Stack http://www.gdargaud.net/Hack/RotAndStack.html
  - DarkTable http://www.darktable.org/
  - RawTherapee http://rawtherapee.com/
  
  (c) E. Farhi, 2018. GPL2.
