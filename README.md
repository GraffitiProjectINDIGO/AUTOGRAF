![Alt text](/images/AUTOGRAF_logo/INDIGO_logoAUTOGRAF.png?raw=true "Optional Title")
## Short Description
AUTOGRAF (**AUT**omated **O**rthorectification of **GRAF**fiti photos) is an open-source python-based Metashape add-on which enables the automated orthorectification of graffiti at a specific site of interest. It employs state-of-the art photogrammetric computer vision techniques to allow highly accurate georeferencing and orthorectification of large numbers of photographs. A paper detailing AUTOGRAF's methodology will soon be submitted to Heritage (an MDPI journal). 

AUTOGRAF is developed as part of the [INDIGO project](https://projectindigo.eu/) (In-ventory and DI-sseminate G-raffiti along the d-O-naukanal) carried out by the [Ludwig Boltzmann Institute for Archaeological Prospection and Virtual Archaeology](https://archpro.lbg.ac.at/) in close collaboration with the [GEO Department of TU Wien University](https://www.geo.tuwien.ac.at/).

## How to set up AUTOGRAF 
Before AUTOGRAF can be used, the following preparatory steps [1-3] need to be performed: 
### 1 - Install Agisoft's Metashape
Agisoft's Metashape version 1.8.3 (earlier version might not be working as AUTOGRAF was designed and tested on 1.8.3) must be installed and the license must be active. More info here: https://www.agisoft.com/downloads/installer/

### 2 - Install external python packages
Some external python packages must installed into Metashape's python environment. These packages are: *numpy*, *matplotlib*, *scikit-image*. To do this the following command must be executed via the command line (note that programme paths might need to be adapted)

#### on Windows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%programfiles%\Agisoft\Metashape Pro\python\python.exe" -m pip install numpy matplotlib scikit-image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### on MacOS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
/MetashapePro.app/Contents/Frameworks/Python.framework/Versions/3.8/bin/python3.8 -m pip install numpy matplotlib scikit-image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#### on Linux 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
./metashape-pro/python/bin/python3.8 -m pip install numpy matplotlib scikit-image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

### 3 - Add AUTOGRAF to Metashape
First you have to download the AUTOGRAF.py script from this repository (./src/AUTOGRAF.py). Then you must add it to the METASHAPE GUI. There are several ways to do this. For example as follows (a or b): 

a) In the Metahsape GUI open *Tools -> Run Script -> select the AUTOGRAF.py*

b) If you use AUTOGRAF regularly, you can also add it permanently to Metashape by copying AUTOGRAF.py to the following directory (this is a Windows example, MacOS and Linux are similar):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 C:/Users/<username>/AppData/Local/Agisoft/Metashape Pro/scripts/
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

Both (a and b) should result in a GUI looking like this:

![Alt text](/images/1.png?raw=true "Optional Title")

## Using AUTOGRAF

AUTOGRAF follows the methodology introduced by Wild et al. (in preparation)[^1]. The general workflow is summarised in the following workflow chart[^1]: 
![Alt text](/images/2.png?raw=true "Optional Title")

### Existing network of oriented images
A prerequisite of the currently implemented version of AUTOGRAF is that a network of oriented cameras exists. This network must be acquired at the studied site prior to using AUTOGRAF as AUTOGRAF uses the existing network to incrementally orient and add new photographs. In Metashape terminology this network is realised through a "main" chunk containing all oriented images. **IMPORTANT:** this "main" chunk must be the first chunk in the list of chunks: 
<p align="center">
<img src="/images/3.PNG?raw=true" alt="Sublime's custom image"/> 
</p>

### How to pass photos to AUTOGRAF
As input AUTOGRAF expects one folder, containing all subfolders that need to be processed. A subfolder contains the images of one graffito: 
<p align="center">
<img src="/images/5.png?raw=true" alt="Sublime's custom image"/> 
</p>

### Run AUTOGRAF
AUTOGRAF is started by clicking the "Run" button in the dropdown menu. AUTOGRAF automatically executes the above methodology and writes the resulting orthophotos in a "results" folder of each subfolder. AUTOGRAF produces 2 main results: 

1. Orthophotos in different raster cell sizes (1cm / 1mm / native)
2. A custom-made georeferencing file (.csv). This file contains the transformation parameters that allow assigning a 3D world coordinate to each pixel. Please note that this file is not standardised (as no standard for this exists yet) 

### Please Note
Currently, AUTOGRAF still has the following limitations: 
- only Metashape version 1.8 is supported
- it only supports rectilinear lenses (e.g. fiesheye lense are NOT supported yet)
- only the Austrian MGI / GK EAST CRS is supported (EPSG: 31256)

All those things can be easily altered in the code for your own use but we will consider implementing above-mentioned features in the future (your feedback is thus highly appreciated!) Please help us improving AUTOGRAF by raising an issue here on GitHub or directly write an E-Mail to bewild@projectindigo.eu

-----------------------------------------------------------------------------------------------------------------------
### How to cite AUTOGRAF
For citing the methodology:
Wild, B., Verhoeven, G., Schlegel, J., Wogrin, S., Wieser, M., Ressl, C., Otepka-Schremmer, J., Pfeifer, N. 2022. AUTOGRAF - AUTomated Orthorectification of GRAFfiti photos. Heritage. Submitted.

if you use the provided code please also cite this:
Wild B. 2022. AUTOGRAF (AUTomated Orthorectification of GRAFfiti photos). Zenodo. https://doi.org/10.5281/zenodo.7049950

Wild, B., Verhoeven, G., Schlegel, J., Wogrin, S., Wieser, M., Ressl, C., Otepka-Schremmer, J., Pfeifer, N. 2022. AUTOGRAF - AUTomated Orthorectification of GRAFfiti photos. Heritage. Submitted.
