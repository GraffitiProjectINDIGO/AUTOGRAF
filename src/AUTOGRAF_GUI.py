# AUTOGRAF - AUTomated Orthorectification of GRAFfiti photos
# ----------------------------------------------------------------------------------------------------------------------
# AUTOGRAF was developed in the framework of the graffiti research project INDIGO (https://projectindigo.eu/).
# INDIGO is funded by the Heritage Science Austria programme of the Austrian Academy of Sciences (ÖAW)
# ----------------------------------------------------------------------------------------------------------------------
# AUTOGRAF allows the fully automated retrieval of graffiti orthophotos from images.
# It is implemented as add-on to Agisoft's Metashape and can be added as toolbox/script to the Metashape GUI. Detailed
# instructions on how to set AUTOGRAF up can be found on the corresponding GitHub repository and information on
# whole methodology can be found in Wild et al. (in preparation for MDPI-Heritage)
# ----------------------------------------------------------------------------------------------------------------------

import os, math, statistics, datetime
from PySide2 import QtGui, QtCore, QtWidgets
import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import ransac, LineModelND
from os.path import isfile, join
from os import listdir
import pickle
from PySide2 import QtGui, QtCore, QtWidgets

global metashapeApp
metashapeApp = Metashape.Application()

# Checking Metashape compatibility #
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #
compatible_major_version = "1.8"
found_major_version = ".".join(metashapeApp.version.split('.')[:2])
if found_major_version != compatible_major_version:
    raise Exception("Incompatible Metashape version: {} != {}".format(found_major_version, compatible_major_version))

# All individual functions #
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #

def parse_world_file(ortho_world_file):
    """
    Helper function for reading world files (.jgw/.tfw/...)
    :param ortho_world_file: path to world file (as exported from Metashape)
    :return: the extracted georeferencing info as list
    """

    with open(ortho_world_file, 'r') as fd:
        lines = fd.read().splitlines()
    geo=[float(x) for x in lines]
    return geo

def estimate_geo_DATA(ortho_rotatoin,ortho_geotransform,point_on_plane):
    """
    Computes the 2D->3D(ortho) transformation parameters
    :param ortho_rotatoin: rotation parameters for 2D-3D(ortho) transformation (3x3)
    :param ortho_geotransform: parameters from world file (as derived from parse_world_file()
    :param point_on_plane: a point within the plane (in 3D world). (e.g.: centroid of sparse PC from plane derivation)
    :return: the 2D-3D (ortho) transformation parmaeters
    """

    # Prepare vector in orthophoto system
    p = np.array(point_on_plane)
    # actually only height of the plane in ortho SYS must be estimated.
    p1_orthoSYS = np.matmul(ortho_rotatoin, p)
    translation = [ortho_geotransform[4], ortho_geotransform[5], p1_orthoSYS[2]]

    return translation, ortho_geotransform[0]

def writeGeoreferenceFile(ortho_R, wld_file_path, point_on_plane, pathToResults, fname):
    """
    Writes a custom georeferencing file which stores the required information for transforming the 2D ortho in the 3D-
    space (from the ortho plane). NOTE: this is not standardised as so far no 3D world files or similar exist
    :param ortho_R: transformation matrix (3x3)
    :param wld_file_path: (path to the Metashape-derived world file)
    :param point_on_plane: a (arbitrary) point (in 3D world-system) within the plane (this code uses the centroid of the sparse point
    cloud.
    :param pathToResults: path to where the file should be stored (i.e. results folder in graffito folder)
    :param fname: filename of the world file ("uniqueGraffitoID".txt)
    :return: -
    """
    wld_file_data = parse_world_file(wld_file_path)
    ortho_T, ortho_m = estimate_geo_DATA(ortho_R, wld_file_data, point_on_plane)
    pathForTxt = os.path.join(pathToResults, fname)
    with open(pathForTxt, 'w') as f:
        f.write(str(ortho_T[0]) + ', ' + str(ortho_T[1]) + ', ' + str(ortho_T[2]))
        f.write('\n' + str(ortho_m))
        f.write('\n' + str(ortho_R[0][0]) + ', ' + str(ortho_R[0][1]) + ', ' + str(ortho_R[0][2]) +
                '\n' + str(ortho_R[1][0]) + ', ' + str(ortho_R[1][1]) + ', ' + str(ortho_R[1][2]) +
                '\n' + str(ortho_R[2][0]) + ', ' + str(ortho_R[2][1]) + ', ' + str(ortho_R[2][2]))
        f.write('\n' + str(point_on_plane[0]) + ', ' + str(point_on_plane[1]) + ', ' + str(point_on_plane[2]))

def getPath():
    """
    Reads the directory path which was set for processing.
    :return: Path to the directory to be processed (containing all graffiti subfolder - each with unique ID/name)
    """
    global pathToGraffitoImages
    pathToGraffitoImages = Metashape.app.getExistingDirectory('Choose directory to be processed')
    return pathToGraffitoImages

def folders_in(path_to_parent):
    """
    Returns all subfolder of the specified directory (i.e. all graffiti folders)
    :param path_to_parent: Parent directory (i.e. folder containing all graffiti folders)
    :return: List of paths to all folders
    """
    for fname in os.listdir(path_to_parent):
        if os.path.isdir(os.path.join(path_to_parent,fname)):
            yield os.path.join(path_to_parent,fname)

def getLastDirectories(subfolderList):
    """
    Returns the names of the folders (not the whole path but only the name!)
    :param subfolderList:
    :return:
    """
    folderGroupNames = []
    for folder in subfolderList:
        folderGroupNames.append(str(folder.split('\\')[1]))
    return folderGroupNames

def getImagePathsFromFolder(folder):
    """
    Extracts the paths to all images in a graffito folder
    :param folder: current graffito folder
    :return: List of paths to graffito-photos
    """
    image_list = [f for f in listdir(folder) if isfile(join(folder, f))]
    photo_list = list()
    global fname
    for photo in image_list:
        if photo.rsplit(".", 1)[1].lower() in ["jpg"]:
            fname = photo
            photo_list.append("/".join([folder, photo]))
    return photo_list

def transform2Dto3DCamera(chunk, point, cam= None):
    """
    Computes 3D world coordinates from 2D pixel coodrinates (for images with known orientation). This is done by
    intersecting the image ray (which goes thorugh the camera center and the defined pixel with the georeferenced 3D
    model). This is primarily needed to perform to define the graffito outline in 3D (needed for selecting only the
    part of interest of the graffito)
    :param chunk: Metashape chunk of the currently processed graffito.
    :param point: the point which should be transformed (2D->3D)
    :param cam: Metashape camera of the camera in which the pixel coordinntes where picked
    :return: The point (3D-wolrd) of the pciked image coordinates
    """
    if cam == None:
        cam = chunk.cameras[-1]
    point_internal = chunk.model.pickPoint(cam.center, cam.unproject(point))
    point3D_referenced = chunk.crs.project(chunk.transform.matrix.mulp(point_internal))
    return point3D_referenced

def intersectWithPlane(chunk, camera, pt_2D, R, C):
    """
    This is the backup function for transform2Dto3DCamera(). It is executed when no 3D model is available for a certain
    image coordinate. This occurs frequently because the graffito outline is often at the very border at of the object
    making it likely that no 3D model can be intersected. This function does not use the 3D model but the RANSAC-derived
    plane which approximates the surface.
    :param chunk: Metashape chunk of the currently processed graffito.
    :param camera: Metashape camera of the camera in which the pixel coordinntes where picked
    :param pt_2D: the point which should be transformed (2D->3D)
    :param R: Rotation matrix (image system - world system; 3x3 matrix)
    :param C: centroid of the plane as derived from the sparse point cloud
    :return: Intersection (3D point) between (camera center - pixel in image x approximating plane of surface)
    """
    C_camera = chunk.crs.project(chunk.transform.matrix.mulp(camera.center))

    img_3d = camera.unproject(pt_2D)
    ray = chunk.crs.project(chunk.transform.matrix.mulp(img_3d)) - chunk.crs.project(
        chunk.transform.matrix.mulp(camera.center))
    b = ray / ray.norm()

    C_plane = Metashape.Vector(C)
    a1 = np.asarray(R).reshape(3, 3)[0, 0:3]
    a2 = np.asarray(R).reshape(3, 3)[1, 0:3]

    A = np.transpose(np.array([a1, a2, -b]))
    param = np.matmul(np.linalg.inv(A), (C_camera - C_plane))

    S1 = a1 * param[0] + a2 * param[1] + C_plane
    # S2 = b * param[2] + C_camera

    return S1


def resizeBoundingBox(chunk):
    """
    automatically resize bounding box based on actual extent of the point cloud. Necessary as this does not happen
    automatically in Metashape when chunks are copied and partly deleted
    :param chunk: chunk which bounding box should be resized.
    :return: -
    """
    region = chunk.region
    T = chunk.transform.matrix

    m = Metashape.Vector([10E+10, 10E+10, 10E+10])
    M = -m

    for point in chunk.point_cloud.points:
        if not point.valid:
            continue
        coord = T * point.coord
        for i in range(3):
            m[i] = min(m[i], coord[i])
            M[i] = max(M[i], coord[i])

    center = (M + m) / 2
    size = M - m
    region.center = T.inv().mulp(center)
    region.size = size * (1 / T.scale())

    region.rot = T.rotation().t()

    chunk.region = region

def filterModelBasedOnConfidence(model, confidenceThreshold=0):
    """
    can be used to filter the model based on the derived confidence values of the mesh (mostly determined by the
    overlap of the images). NOTE: currently not used because it leads to an increased failure rate of 3D graffito
    outline computation
    :param model: 3D mesh
    :param confidenceThreshold: int (see Metashape documentation)
    :return: -
    """
    for face in model.faces:
        for i in face.vertices:
            if model.vertices[i].confidence <= confidenceThreshold:
                face.selected = True
                continue
    model.removeSelection()
    #model.closeHoles() # can be activated to close potential holes.
    for face in model.faces:
        face.selected = False

def removeIsoltatedComponents(model, threshold=0.33):
    """
    Removes component of the 3D mesh which are not connected to the "main" model. These are usually branches and
    other vegetation/infrastructure.
    :param model: 3D mesh
    :param threshold: percentage of number of vertices the component can have maximally to not be removed (i.e. if 33%
    every component which has less than 33% of the total number of vertices will be removed).
    :return: -
    """
    faceCount = len(model.faces)
    model.removeComponents(int(faceCount*threshold))

def ransacCleanMesh(model):
    """
    To still be fully implemented!
    This function should (in the future) allow a more robust mesh cleaning (not only based on aboves
    remove IsolatedComponents() but also taking into account the distance of the vertices to the RANSAC derived
    plane.
    :param model: 3D mesh
    :return: - NOT YET FULLY IMPLEMENTED
    """
    X,Y = [], []
    # vertices = model.vertices
    for vertex in model.vertices:
        X.append(vertex.coord[0])
        Y.append(vertex.coord[1])
    data = np.column_stack([X, Y])
    modelRANSAC = LineModelND()
    modelRANSAC.estimate(data)
    model_robust, inliers = ransac(data, LineModelND, min_samples=2, residual_threshold=0.5, max_trials=1000)
    outliers = inliers == False
    P1 = modelRANSAC.params[0] - 100 * modelRANSAC.params[1]
    P2 = modelRANSAC.params[0] + 100 * modelRANSAC.params[1]

    vertices = model.vertices
    faces = model.faces

    for face in faces:
        face_vertices = [vertices[v] for v in face.vertices]

        fit_vert = 0
        for vertex in face_vertices:
            P3 = vertex.coord[0:2]
            if (np.linalg.norm(np.cross(P2-P1, P1-P3))/np.linalg.norm(P2-P1)) > 0.5:
                face.selected = True

    model.removeSelection()

    # generate coordinates of estimated models
    line_x = np.arange(np.quantile(X, 0.1)-1, np.quantile(X, 0.9)+1)

    return inliers

def selectCamerasBasedonRadius(inputCamCoord, existingCams, existingCamChunk, radius=30):
    """
    Selects cameras that are within a certain radius to a point. This helps to remove the majority of cameras and
    significantly boost processing times.
    :param inputCamCoord: camera coordinates of all oriented cameras for one graffito. These are used to compute
    a ficitive "median" camera pose.
    :param existingCams: all cameras in the block
    :param existingCamChunk: chunk which contains all cameras (the "major chunk" containing the whole image network)
    :param radius: radius (in metres) beyond which cameras are deactivated and thus ignored for the incremental SfM
    :return: -
    """
    crsWGS84 = Metashape.CoordinateSystem("EPSG::4326")
    crsMGI = Metashape.CoordinateSystem("EPSG::31256")

    inputCamCoordMGI = []
    for coord in inputCamCoord:
        inputCamCoordMGI.append(Metashape.CoordinateSystem.transform(Metashape.Vector(coord), crsWGS84, crsMGI))

    inputCoordMeanMGI = np.median(np.asarray(inputCamCoordMGI), axis=0)

    for cam in existingCams:
        if cam.transform:
            try:
                camCenter = existingCamChunk.crs.project(existingCamChunk.transform.matrix.mulp(cam.center))
                dist = math.dist(inputCoordMeanMGI[0:2], camCenter[0:2])

                if dist < searchRadius:
                    cam.enabled = True
            except:
                print(cam)


def initialBundleblock():
    """
    It computes an intial local bundleBlock of the currently processed graffito folder. The main results form this are:
        - the inner orientation of the camera (iOr)
        - the exterior orientations (exOr)
        - the sparse cloud

    This function also reads the GNSS information from a potentially mounted GNSS receiver. This info is written into
    a pickle file and stored in the folder where the graffito-block images are located.

    Most important parameters:
        - key point limit: 25000
        - tie point limit: 3000
    :return:
    """
    subfolderList = list(folders_in(pathToGraffitoImages)) # get all subfolders from the selected directory
    folderNames = getLastDirectories(subfolderList)

    for currentFolder, folderName in zip(subfolderList, folderNames):
        photo_list = getImagePathsFromFolder(currentFolder)

        chunk = doc.addChunk()
        chunk.label = 'Initial SfM - ' + folderName # name chunk like last file name in selected directory
        print(photo_list)
        chunk.addPhotos(photo_list) # add all images that are in the selected directory

        # If the image quality differs (e.g. due to automated acquisition of images on a moving platform) the
        # following three lines of code can be used to check the image quality before processing the whole block
        # boolean = Metashape.app.getBool('Do you want to check the image quality before starting the process?')
        # if boolean:
        #    checkImageQuality(chunk)
        camEstLoc = []
        for cam in chunk.cameras:
            camEstLoc.append(np.asarray(cam.reference.location))

        estCamLocFn = '/'.join(chunk.cameras[-1].photo.path.split('/')[:-1]) + '/' + chunk.label[
                                                                                          14:] + '_GNSS_coord.pickle'

        with open(estCamLocFn, 'wb') as handle:
            pickle.dump(camEstLoc, handle, protocol=pickle.HIGHEST_PROTOCOL)

        chunk.matchPhotos(downscale=dsImageMatching, generic_preselection=True, reference_preselection=True, keep_keypoints=True, tiepoint_limit=tiePoints, keypoint_limit=keyPoints) # perform initial bundleblock
        chunk.alignCameras(reset_alignment=False) # perform initial bundleblock

        # save initial IOR estimation in seperate file for potential later use
        sensor = chunk.sensors[0]
        sensor.calibration.save('/'.join(chunk.cameras[-1].photo.path.split('/')[:-1]) + '/' + 'iOr.xml')

def initialBundleBlockChecks():
    '''
    Here, the results from the initial bundleblock are analysed and validated. Mainly this functions does the following:
        - checking if all images were oriented. (Filenames of) images which were not oriented are written into .txt
          files. This (potential) list of files is used in later steps to deactivate images which should be ignored
          in the further processing steps. This decreases processing times and increases reliability.
        - computed various error metrics from the initial bundleblock. Most importantly the RMSE_img is derived. The
          RMSE_img is a per image representation of the reprojection error. The RMSE_img is written into a .txt file and
          stored along with the graffito-block images. No action is automatically executed based on these metric
          because threshold-based removal of images might result in holes in the mesh and thus failure of the
          subsequent workflow. However, the RSME_img can be used to manually intervene at a later stage.
          (-> see getPointErrors() for more details)

    :return: The function writes 4 files:
                - "uniqueID"_failed.pickle: stores filenames of the images which could not be aligned during the intial
                  SfM.
                - "uniqueID"_inaccurate.pickle: stores filename of images for which the RMSE_img exceeds a certain
                  threshold. This threshold can be set in getPointErrors(). Default is 1.5 (pixels)
                - "uniqueID"_bad.pickle: union of the two above.
    '''

    subfolderList = list(folders_in(pathToGraffitoImages))
    folderGroupNames = getLastDirectories(subfolderList)

    for chunk in doc.chunks[1:]:
        if chunk.label[14:] in folderGroupNames:
            failedCameras, goodCameras, badCameras = [], [], []
            print('checking for: ' + chunk.label)
            for camera in chunk.cameras:
                if camera.transform:
                    goodCameras.append(camera.photo.path)
                else:
                    print("\n !!!!!!!!!!!!! " + camera.label + " not aligned successfully !!!!!!!!!!! \n")
                    failedCameras.append(camera.photo.path)

            inaccurateCameras = getPointErrors(chunk) # calculates error statistics
            for inaccurateCamera in inaccurateCameras:
                badCameras.append(inaccurateCamera)
            for failedCamera in failedCameras:
                badCameras.append(failedCamera)
            if len(badCameras) == 0:
                print('All cameras are OK')
            else:
                print('the following cameras are problematic and will be removed: ')
                for failedCamera in failedCameras:
                    print(failedCamera)

            for camera in chunk.cameras:
                if camera.photo.path not in badCameras:
                    goodCameras.append(camera)

            pickle_fname_bad = '/'.join(chunk.cameras[-1].photo.path.split('/')[:-1]) + '/' + chunk.label[14:] + '_bad.pickle'
            pickle_fname_failed = '/'.join(chunk.cameras[-1].photo.path.split('/')[:-1]) + '/' + chunk.label[14:] +  '_failed.pickle'
            pickle_fname_inaccurate = '/'.join(chunk.cameras[-1].photo.path.split('/')[:-1]) + '/' + chunk.label[14:] + '_inaccurate.pickle'
            with open(pickle_fname_bad, 'wb') as handle:
                pickle.dump(badCameras, handle, protocol=pickle.HIGHEST_PROTOCOL)
            with open(pickle_fname_failed, 'wb') as handle:
                pickle.dump(failedCameras, handle, protocol=pickle.HIGHEST_PROTOCOL)
            with open(pickle_fname_inaccurate, 'wb') as handle:
                pickle.dump(inaccurateCameras, handle, protocol=pickle.HIGHEST_PROTOCOL)

            print(badCameras)
            print(inaccurateCameras)
            print(failedCameras)

    print('All checks performed. Please check the Metashape GUI and selected graffiti folders for the results.')



def getPointErrors(chunk):
    '''
    Copmutes various error metrics from the initial SfM. Main result is the the RMSE_img (as described in
    initialBundleBlockChecks())
    :param chunk: the currently processed chunk (i.e. the current graffito block)
    :return: list of cameras which exceeds a certain threshold for the RMSE_img. Also the error metrics are stored
    seperately in the folder of the respective currently processed graffito block. Lastly, additional quality
    is printed to the console.
    '''

    # Compute the re-projection error
    point_cloud = chunk.point_cloud
    points = point_cloud.points
    error, tie_points = [], []

    badCameras = []

    camera_metrics = {}

    for camera in chunk.cameras:
        if camera.transform:
            error_temp = []
            point_index = 0
            photo_num = 0
            for proj in chunk.point_cloud.projections[camera]:
                track_id = proj.track_id
                while point_index < len(points) and points[point_index].track_id < track_id:
                    point_index += 1
                if point_index < len(points) and points[point_index].track_id == track_id:
                    if not points[point_index].valid:
                        continue

                    dist = camera.error(points[point_index].coord, proj.coord).norm() ** 2
                    error.append(dist)
                    error_temp.append(dist)
                    photo_num += 1
                    tie_points.append(photo_num)
            camera_metrics[camera.label] = round(math.sqrt(sum(error_temp) / len(error_temp)), 5)#np.mean(np.asarray(error_temp))#

    with open(chunk.cameras[-1].photo.path[:-4] + "_repr_Error" + ".txt", "wt") as f:
        print(camera_metrics, file=f)

    for camera in chunk.cameras:
        if camera.transform:
            print('------------------------------------------------')
            print('RMS reprojection error for ' + camera.label + ' = ' + str(camera_metrics[camera.label]))
            if camera_metrics[camera.label] > 1.5:
                print(camera.label + ' is considered bad and will be removed')
                badCameras.append(camera.photo.path)
                #chunk.remove(camera)
            print('------------------------------------------------')
    f.close()

    reprojection_rmse = round(math.sqrt(sum(error) / len(error)), 2)
    reprojection_max = round(max(error) , 2)
    reprojection_std = round(statistics.stdev(error), 2)
    tie_points_per_image = round(sum(tie_points) / len(tie_points), 0)


    print("Average RMS reprojection error: " + str(reprojection_rmse))
    #print("Average tie point residual error: " + str(np.mean(np.asarray(error))))
    #print("Maxtie point residual error: " + str(np.max(np.asarray(error))))
    print("Max tie RMS reprojection error: " + str(reprojection_max))
    print("Standard deviation for tie point residual error: " + str(reprojection_std))
    #print("Standard deviation for tie point residual error: " + str(np.std(np.asarray(error))))
    print("Average number of tie points per image: " + str(tie_points_per_image))

    #fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    #plt.title('Histogram of RMS reprojection error - all tie points')
    #axs.hist(error, bins=np.linspace(0,1,200))
    #plt.axvline(np.mean(np.asarray(error)), color='k', linestyle='dashed', linewidth=1)
    #plt.show()
    #input("Press Enter to continue...")
    #plt.close()

    return badCameras

def checkImageQuality(chunk):
    '''
    checks the image quality for each photo in a chunk (using the checkImage routine implemented in Metashape)
    :param chunk: chunk to be processed
    :return:
    '''
    lowQualityImages = []
    for image in chunk.cameras:
        chunk.analyzePhotos(image)
        if float(image.meta['Image/Quality']) < 0.6:
            lowQualityImages.append(image)

    if len(lowQualityImages) != 0:
        print('The following images have a image quality of below 0.6:')
        for image in lowQualityImages:
            print(image.label)
        boolean = Metashape.app.getBool(str(len(lowQualityImages)) + ' images have a low quality. Remove those images (y - recommended) '
                                 'or keep (n) those images')
        if boolean:
            for image in lowQualityImages:
                chunk.remove(image)
                print(image.label + ' removed')


def incrementalSfM(downScalingImageMatching):
    '''
    This function adds all images, that should be incrementally oriented to the "main" chunk. (NOTE: the "main" chunk
    is always defined as the first chunk in the chunk list!). Once all photos are added they are incrementally
    oriented. For speeding up the computations, available GNSS data is used to limit the search radius. Note that
    all images are oriented in one go, so this process can take several hours depending on the number of images (see
    Wild et al, 2022 for details regarding processing times).
    :param downScalingImageMatching: int, determines the downscaling that is applied to the images before the image
    matching procedure starts. The higher the number, the lower is the images' resolution resulting in faster
    prcessing but less accurate results. possible values: [0.25, 1 (no downscaling is applied), 4, 16, 64]
    '''
    chunk = doc.chunks[0]
    for cam in chunk.cameras:
        cam.enabled= False

    # set up the camera groups in the main chunk for subalignment
    cameraGroupsExisting, cameraGroupLabelsOriginal  = chunk.camera_groups, []
    for camGroup in cameraGroupsExisting:
        cameraGroupLabelsOriginal.append(camGroup.label)

    cameraGroupsAdded = []
    subfolderList = list(folders_in(pathToGraffitoImages))
    existingCameras = chunk.cameras

    sensorNames = []
    for chunkSensors in chunk.sensors:
        sensorNames.append(chunkSensors.label)

    for folder in subfolderList:
        #parentFolder = str(folder.split('/')[-1])
        folderGroupName = str(folder.split('\\')[1])
        cameraGroupsAdded.append(folderGroupName)
        if folderGroupName in cameraGroupLabelsOriginal:
            continue
        folderGroup = chunk.addCameraGroup()
        folderGroup.label = folderGroupName
        photo_list = getImagePathsFromFolder(folder)

        chunk.addPhotos(photo_list)
        newCams = []
        sensor = chunk.addSensor()
        sensor.label = folderGroupName

        for camera in chunk.cameras:
            if camera.group is None:
                #sensor.calibration.load('/'.join(camera.photo.path.split('/')[:-1]) + '/' + 'iOr.xml')
                #print('/'.join(camera.photo.path.split('/')[:-1]) + '/' + 'iOr.xml')
                #sensor = chunk.sensors[0]
                sensor.type = Metashape.Sensor.Type.Frame
                sensor.calibration = camera.sensor.calibration
                sensor.width = camera.sensor.width
                sensor.height = camera.sensor.height
                sensor.focal_length = camera.sensor.focal_length
                sensor.pixel_height = camera.sensor.pixel_height
                sensor.pixel_width = camera.sensor.pixel_width

                camera.sensor = sensor

                tempCam = camera
                newCams.append(camera)
                camera.group = folderGroup

        estCamLocFn = '/'.join(tempCam.photo.path.split('/')[:-1]) + '/' + folderGroupName +'_GNSS_coord.pickle'
        with open(estCamLocFn, "rb") as input_file:
            newCamCoord = pickle.load(input_file)

        print(newCamCoord)

        selectCamerasBasedonRadius(newCamCoord, existingCameras, chunk)

    chunk.matchPhotos(downscale=dsImageMatching, generic_preselection=True, reference_preselection=True, keypoint_limit=keyPoints, tiepoint_limit= tiePoints, keep_keypoints=True)
    chunk.alignCameras(reset_alignment=False)

def prepareOrtho():
    '''
    routine to prepare and clean up the folder/chunk sturcture in Metashape. Mainly, this function makes a copy of
    the "main" chunk and creates an individual chunk for each inputted graffito-block. The name of the respective chunk
    is its unique ID.
    :return:
    '''
    mainChunk = doc.chunks[0]
    tempChunk = mainChunk.copy(items=[], keypoints=False)
    tempChunk.label = "Temporary Chunk"

    subfolderList = list(folders_in(pathToGraffitoImages))
    folderGroupNames = getLastDirectories(subfolderList)

    for cameraGroup in tempChunk.camera_groups:
        print(cameraGroup)
        if cameraGroup.label not in folderGroupNames:
            tempChunk.remove(cameraGroup)

    for cameraGroup, i in zip(tempChunk.camera_groups, range(0, len(tempChunk.camera_groups))):
        individualChunk = tempChunk.copy(items=[], keypoints=False)
        individualChunk.label = cameraGroup.label
        for camGroup in individualChunk.camera_groups:
            if camGroup.label != folderGroupNames[i]:
                individualChunk.remove(camGroup)

        badPhotosPath = os.path.join(subfolderList[i], folderGroupNames[i]+'_failed.pickle')
        with open(badPhotosPath, "rb") as input_file:
            badPhotosList = pickle.load(input_file)
        for camera in individualChunk.cameras:
            if camera.photo.path in badPhotosList:
                camera.enabled = False

    doc.remove(tempChunk)

def createOrthos():
    '''
    This function executes all last steps towards the orthocreation. Those are:
        - creation of depth maps
        - 3D model derivation
        - projection plane computation
        - orthophoto computation
        - segmentation of the orthophoto (if segmentation data is available)
    :return:
    '''
    ### prepare orthorectification
    subfolderList = list(folders_in(pathToGraffitoImages))
    folderGroupNames = getLastDirectories(subfolderList)
    i = 0

    for chunk in doc.chunks:
        try:
            if chunk.label not in folderGroupNames:
                continue
            else:
                counter = 0
                for camera in chunk.cameras:
                    if camera.enabled == False:
                        counter = counter + 1
                if counter == len(chunk.cameras):
                   #'For chunk:' + chunk.label + ' all cameras were diabled because they did not pass the'
                    #                                 'accuracy checks. Continuing with next chunk.')
                    i = i + 1
                    continue
                #sensor = chunk.sensors[0]
                #sensor.calibration.load('/'.join(chunk.cameras[-1].photo.path.split('/')[:-1]) + '/' + 'iOr.xml')

                resizeBoundingBox(chunk)
                R, C = getProjectionPlane(chunk)
                orthoproj = Metashape.OrthoProjection()
                orthoproj.crs = chunk.crs
                orthoproj.matrix = Metashape.Matrix.Rotation(R)

                chunk.buildDepthMaps(downscale=dsDenseMatching)

                #### If the depth map creation fails because of GPU-compability issues (which occur with older GPUs
                # sometimes) the following code block can be used as workaround (the buildDepthMaps() should be ingored
                # in this case

                #task = Metashape.Tasks.BuildDepthMaps()
                #task.downscale = downScalingDenseMatching
                #task.filter_mode = Metashape.MildFiltering
                #task["pm_enable"] = False
                #task.apply(chunk)

                chunk.buildModel(source_data=Metashape.DataSource.DepthMapsData, surface_type=Metashape.Arbitrary,
                                 interpolation=Metashape.EnabledInterpolation)

                removeIsoltatedComponents(chunk.model)
                # filterModelBasedOnConfidence(chunk.model, confidenceThreshold=4) # might be implemented in the future

                chunk.buildOrthomosaic(surface_data=Metashape.ModelData, projection=orthoproj, cull_faces=True)
                chunk.orthomosaic.label = chunk.label
                resolution = str(chunk.orthomosaic.resolution * 1000)[:4]
                fileName = chunk.label + '_ortho_' + resolution + 'mm.tif'
                fileName1mm = chunk.label + '_ortho_1mm.tif'
                fileNamePNG = chunk.label + '_ortho_1cm.png'

                # segmentation
                pathToSegment = os.path.join(subfolderList[i], "img_segmented.txt")
                pathToResults = os.path.join(subfolderList[i], "results")
                if not os.path.exists(pathToResults):
                    os.mkdir(pathToResults)
                pathForOrtho = os.path.join(pathToResults, fileName)
                pathForOrtho1mm = os.path.join(pathToResults, fileName1mm)
                pathForOrthoPNG = os.path.join(pathToResults, fileNamePNG)
                shpFname = chunk.label + '_polygon.shp'
                pathForShape = os.path.join(pathToResults, shpFname)
                if os.path.isfile(pathToSegment):
                    try:
                        original_array = np.loadtxt(pathToSegment)  # .reshape(4, 2)
                        reshaped_array = original_array.reshape(int(len(original_array) / 2), 2)

                        pointsSegmented = []
                        for camera in chunk.cameras:
                            if "segment" in camera.label:
                                segmentCam = camera
                        for row in reshaped_array:
                            point2D = Metashape.Vector(row)
                            try:
                                point3D = transform2Dto3DCamera(chunk, point2D, segmentCam)
                            except:
                                point3D = intersectWithPlane(chunk, segmentCam, point2D, R, C)
                                print(point3D)
                            pointsSegmented.append(point3D)

                        # create shape
                        if len(pointsSegmented) > 2:
                            if not chunk.shapes:
                                chunk.shapes = Metashape.Shapes()
                            chunk.shapes.remove(chunk.shapes)
                            chunk.shapes.crs = chunk.crs

                            shape = chunk.shapes.addShape()
                            shape.label = segmentCam.label
                            shape.attributes["Photo"] = camera.label
                            # shape.group = segmented
                            shape.geometry = Metashape.Geometry.Polygon(pointsSegmented)
                            shape.boundary_type = Metashape.Shape.OuterBoundary
                            chunk.exportShapes(path=pathForShape, save_polygons=True)
                        chunk.exportRaster(path=pathForOrtho, source_data=Metashape.OrthomosaicData, save_world=True,  clip_to_boundary=True)
                        writeGeoreferenceFile(np.asarray(R).reshape(3,3), os.path.join(pathToResults, chunk.label + '_ortho_' + resolution + 'mm.tfw'), C, pathToResults, chunk.label + '_ortho_' + resolution + 'mm_georef.txt')
                        chunk.exportRaster(path=pathForOrtho1mm, source_data=Metashape.OrthomosaicData, resolution_x= 0.001, save_world=True, resolution_y = 0.001, clip_to_boundary=True)
                        writeGeoreferenceFile(np.asarray(R).reshape(3,3), os.path.join(pathToResults, chunk.label + '_ortho_1mm.tfw'), C, pathToResults, chunk.label + '_ortho_1mm_georef.txt')
                        chunk.exportRaster(path=pathForOrthoPNG, source_data=Metashape.OrthomosaicData, resolution_x= 0.01, save_world=True, resolution_y = 0.01, clip_to_boundary=True)
                        writeGeoreferenceFile(np.asarray(R).reshape(3,3), os.path.join(pathToResults, chunk.label + '_ortho_1cm.pgw'), C, pathToResults, chunk.label + '_ortho_1cm_georef.txt')
                    except Exception as e:
                        print(e)
                        chunk.exportRaster(path=pathForOrtho, source_data=Metashape.OrthomosaicData,save_world=True,  clip_to_boundary=False)
                        writeGeoreferenceFile(np.asarray(R).reshape(3,3), os.path.join(pathToResults, chunk.label + '_ortho_' + resolution + 'mm.tfw'), C, pathToResults, chunk.label + '_ortho_' + resolution + 'mm_georef.txt')
                        chunk.exportRaster(path=pathForOrtho1mm, source_data=Metashape.OrthomosaicData, clip_to_boundary=False, save_world=True, resolution_x= 0.001, resolution_y = 0.001)
                        writeGeoreferenceFile(np.asarray(R).reshape(3,3), os.path.join(pathToResults, chunk.label + '_ortho_1mm.tfw'), C, pathToResults, chunk.label + '_ortho_1mm_georef.txt')
                        chunk.exportRaster(path=pathForOrthoPNG, source_data=Metashape.OrthomosaicData, save_world=True, resolution_x=0.01, resolution_y=0.01, clip_to_boundary=False)
                        writeGeoreferenceFile(np.asarray(R).reshape(3,3), os.path.join(pathToResults, chunk.label + '_ortho_1cm.pgw'), C, pathToResults, chunk.label + '_ortho_1cm_georef.txt')

                else:
                    chunk.exportRaster(path=pathForOrtho, source_data=Metashape.OrthomosaicData, save_world=True)
                    writeGeoreferenceFile(np.asarray(R).reshape(3,3), os.path.join(pathToResults, chunk.label + '_ortho_' + resolution + 'mm.tfw'), C, pathToResults, chunk.label + '_ortho_' + resolution + 'mm_georef.csv')
                    chunk.exportRaster(path=pathForOrtho1mm, source_data=Metashape.OrthomosaicData, save_world=True, resolution_x= 0.001, resolution_y = 0.001)
                    writeGeoreferenceFile(np.asarray(R).reshape(3,3), os.path.join(pathToResults, chunk.label + '_ortho_1mm.tfw'), C, pathToResults, chunk.label + '_ortho_1mm_georef.txt')
                    chunk.exportRaster(path=pathForOrthoPNG, source_data=Metashape.OrthomosaicData, save_world=True, resolution_x=0.01, resolution_y=0.01, clip_to_boundary=False)
                    writeGeoreferenceFile(np.asarray(R).reshape(3,3), os.path.join(pathToResults, chunk.label + '_ortho_1cm.pgw'), C,pathToResults, chunk.label + '_ortho_1cm_georef.txt')
                i = i + 1
        except Exception as e:
            i = i+1
            print('!!!!!! FAILED !!!!!!')
            print(e)
            continue



def getProjectionPlane(chunk):
    '''
    computes projection plane onto which the orthophoto is orthogonally projected (via previosuly intersecting the 3D
    model). Part of this function is a RANSAC algorithm (ransacINDIGO()) which runs a RANSAC classifcation on the
    sparse point cloud and sperate inliers and outliers (i.e. points assumed to be part of the graffiti-covered
    surface and points which are not part of the surface). Based on this classification, the plane parameters are
    derived from the (classified) sparse point cloud of the scene via solving an eigenvectors problem. The 3 normalised
    eigenvectors are the main results of this function. They are summarised in a matrix (R)
    :param chunk: chunk to be processed
    :return: * R (3x3) matrix containing the plane defining vectors. The order is the following:
                    - first row: ~ horizontal
                    - second row: ~ vertical
                    - third row: normal to the two vectors above
            * C (1x3) vector containing the centroid of the (inlier) point cloud)

    '''
    for camera in chunk.cameras:
        if camera.transform:
            continue
        else:
            print(camera.label + " not aligned successfully")
            chunk.remove(camera)

    if len(chunk.cameras) == 0:
        return 999
    crs = chunk.crs
    T = chunk.transform.matrix
    points = chunk.point_cloud.points
    X, Y, Z = [], [], []

    camCoord = []
    for cam in chunk.cameras:
        camCoord.append(chunk.crs.project(chunk.transform.matrix.mulp(cam.center)))
    camCoordMean = np.mean(np.asarray(camCoord), axis=0)

    for i in range(len(points)):
        pt = Metashape.Vector((points[i].coord[0], points[i].coord[1], points[i].coord[2]))
        point_temp = crs.project(T.mulp(pt))
        X.append(point_temp[0])
        Y.append(point_temp[1])
        Z.append(point_temp[2])

    inlier = ransacINDIGO(X,Y,Z, np.asarray(camCoord), chunk)

    X, Y, Z = np.asarray(X)[inlier], np.asarray(Y)[inlier], np.asarray(Z)[inlier]
    C = np.array([np.mean(X), np.mean(Y), np.mean(Z)])
    coords = np.array([X, Y, Z]).transpose()

    Q = coords - C
    Qt = Q.transpose()
    [V, D] = np.linalg.eig(np.matmul(Qt, Q))
    idx = np.argsort(V)

    outpoint2 = camCoordMean

    z_v_def = np.array([0, 0, 1])
    a1 = D[:, idx[0]]
    a1 = a1 / np.linalg.norm(a1)


    signum = np.sign((C - outpoint2) * a1)
    a1_n = a1*signum
    horizontal = np.cross(a1_n, z_v_def) / np.linalg.norm(np.cross(a1_n, z_v_def))
    vertical = np.cross(horizontal, a1_n) / np.linalg.norm(np.cross(horizontal, a1_n))
    normal = a1_n
    R = Metashape.Matrix([horizontal, vertical, -normal])

    return R, C

def ransacINDIGO(X,Y,Z, camCoords, chunk):
    '''
    Takes a point cloud and runs a RANSAC calssifation to seperate outliers (i.e. points not belonging to the
    reference plane) from inliers (i.e. points not belonging to the reference plane)
    :param X: X-coordiantes of the PC
    :param Y: Y-coordiantes of the PC
    :param Z: Z-coordiantes of the PC
    :param camCoords: coordinates of the camera for the respective graffito block. Needed to define orientation of the
    plane
    :param chunk: chunk to be processed
    :return: list on (3x1) vectors -> the inlier points
    '''

    data = np.column_stack([X, Y])
    data_all = np.column_stack([X, Y, Z])
    model = LineModelND()
    model.estimate(data)

    # robustly detect plane inlier data with RANSAC algorithm
    model_robust, inliers = ransac(data, LineModelND, min_samples=2, residual_threshold=maxPlaneResid, max_trials=1000)
    outliers = inliers == False

    # generate coordinates of estimated models
    line_x = np.arange(np.quantile(X, 0.1)-1, np.quantile(X, 0.9)+1)
    line_y = model.predict_y(line_x)
    line_y_robust = model_robust.predict_y(line_x)

    ## camera positions
    camPosX = camCoords[:,0]
    camPosY = camCoords[:,1]
    camCoordMean = np.mean(camCoords, axis=0)

    # Plotting routine for showing the ransac results (stored in graffito-folder)
    plt.close('all')
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111)
    ax.plot(data[inliers, 0], data[inliers, 1], '.b', alpha=0.6, label='Inlier')
    ax.plot(data[outliers, 0], data[outliers, 1], '.r', alpha=0.6, label='Outlier')
    ax.plot(camPosX, camPosY, 'xg', label='Camera positions')
    ax.plot(line_x, line_y, '-k', label='Linear Regression')
    ax.plot(line_x, line_y_robust, '-b', label='RANSAC')
    ax.legend(loc='lower left')
    ax.set_xlabel('Rechtswert (X)')
    ax.set_ylabel('Hochwert (Y)')
    ax.set_aspect('equal', 'box')
    ax.set_xlim([np.quantile(X, 0.1) - 5, np.quantile(X, 0.9) + 5])
    ax.set_ylim([np.quantile(Y, 0.1) - 5, np.quantile(Y, 0.9) + 5])
    fig_fname = os.path.dirname(chunk.cameras[-1].photo.path) + '/RANSAC_result.png'
    plt.savefig(fig_fname)

    X_inlier, Y_inlier, Z_inlier = data_all[inliers, 0], data_all[inliers, 1], data_all[inliers, 2]
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(data_all[outliers, 0], data_all[outliers, 1], data_all[outliers, 2], s=5,  alpha=0.5, c='r')
    ax.scatter(X_inlier, Y_inlier, Z_inlier, s=15, c='g')
    ax.set_xlim([np.quantile(X, 0.1)-1, np.quantile(X, 0.9)+1])
    ax.set_ylim([np.quantile(Y, 0.1)-1, np.quantile(Y, 0.9)+1])
    ax.set_zlim([np.quantile(Z, 0.1)-1, np.quantile(Z, 0.9)+1])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    #plt.show()
    plt.close('all')

    return inliers


def dummyFunc():
    return None

def about():
    metashapeApp.messageBox("AUTOGRAF - AUTomated Orthorectification of GRAFfiti photos. \n"
                   "AUTOGRAF was developed in the framework of the graffiti research project INDIGO")

def runAll():
    '''
    Function that executes AUTOGRAF's core functions and passes the parameters to the script (downscaling for image
    matching and depth map generation)
    '''


    print('!!!!  Start: Initial Bundleblock')
    print(datetime.datetime.now())
    initialBundleblock()

    print('!!!!  Initial Bundleblock Checks: ')
    print(datetime.datetime.now())
    initialBundleBlockChecks()

    print('!!!!  Incremental SfM: ')
    print(datetime.datetime.now())
    incrementalSfM(dsImageMatching)

    print('!!!!  Ortho Preparation: ')
    print(datetime.datetime.now())
    prepareOrtho()

    print('!!!!  Ortho Creation: ')
    print(datetime.datetime.now())
    createOrthos()

    doc.save()

    print('!!!!  Finished: ')
    print(datetime.datetime.now())



# Add new items to menu #
# ^^^^^^^^^^^^^^^^^^^^^ #
if __name__ == '__main__':
    '''
    Adds the toolbox to the Metashape GUI
    '''
    plt.close('all')
    global doc
    doc = Metashape.app.document

    class ProgDlg(QtWidgets.QDialog):

        def __init__(self, parent):
            QtWidgets.QDialog.__init__(self, parent)
            self.setWindowTitle("Custom Processing")

            self.btnQuit = QtWidgets.QPushButton("&Exit")
            self.btnP1 = QtWidgets.QPushButton("&Process!")
            self.pBar = QtWidgets.QProgressBar()
            self.pBar.setTextVisible(False)

            self.applyTxt = QtWidgets.QLabel("Accuracy")
            self.radioBtn_high = QtWidgets.QRadioButton("High (slowest)")
            self.radioBtn_medium = QtWidgets.QRadioButton("Medium")
            self.radioBtn_slow= QtWidgets.QRadioButton("Low (fastest)")
            self.radioBtn_high.setChecked(True)
            self.radioBtn_medium.setChecked(False)
            self.radioBtn_slow.setChecked(False)

            # Search radius
            self.intTxt_tiePoints = QtWidgets.QLabel("Tie point limit:")
            self.intEdt_tiePoints = QtWidgets.QLineEdit()
            self.onlyInt = QtGui.QIntValidator()
            self.intEdt_tiePoints.setPlaceholderText("3000")
            self.intEdt_tiePoints.setValidator(self.onlyInt)

            # Search radius
            self.intTxt_keyPoints = QtWidgets.QLabel("Key point limit:")
            self.intEdt_keyPoints = QtWidgets.QLineEdit()
            self.onlyInt = QtGui.QIntValidator()
            self.intEdt_keyPoints.setPlaceholderText("25000")
            self.intEdt_keyPoints.setValidator(self.onlyInt)

            # Search radius
            self.intTxt_searchRadius = QtWidgets.QLabel("Search radius [m]:")
            self.intEdt_searchRadius = QtWidgets.QLineEdit()
            self.onlyInt = QtGui.QIntValidator()
            self.intEdt_searchRadius.setPlaceholderText("30 m")
            self.intEdt_searchRadius.setValidator(self.onlyInt)

            # max plane residuals
            self.intTxt_maxPlaneResid  = QtWidgets.QLabel("Max plane residuals [cm]:")
            self.intEdt_maxPlaneResid = QtWidgets.QLineEdit()
            self.onlyInt = QtGui.QIntValidator()
            self.intEdt_maxPlaneResid.setPlaceholderText("10 cm")
            self.intEdt_maxPlaneResid.setValidator(self.onlyInt)


            layout = QtWidgets.QGridLayout()
            # layout.setSpacing(5)
            layout.addWidget(self.applyTxt, 0, 0)
            layout.addWidget(self.radioBtn_high, 1, 0)
            layout.addWidget(self.radioBtn_medium, 2, 0)
            layout.addWidget(self.radioBtn_slow, 3, 0)
            layout.addWidget(self.intTxt_searchRadius, 0, 3)
            layout.addWidget(self.intEdt_searchRadius, 0, 2)
            layout.addWidget(self.intTxt_maxPlaneResid, 1, 3)
            layout.addWidget(self.intEdt_maxPlaneResid, 1, 2)
            layout.addWidget(self.intTxt_tiePoints, 2, 3)
            layout.addWidget(self.intEdt_tiePoints, 2, 2)
            layout.addWidget(self.intTxt_keyPoints, 3, 3)
            layout.addWidget(self.intEdt_keyPoints, 3, 2)
            layout.addWidget(self.pBar, 4, 0)
            layout.addWidget(self.btnP1, 4, 1)
            layout.addWidget(self.btnQuit, 4, 2)
            self.setLayout(layout)

            proc = lambda: self.process()

            QtCore.QObject.connect(self.btnP1, QtCore.SIGNAL("clicked()"), proc)
            QtCore.QObject.connect(self.btnQuit, QtCore.SIGNAL("clicked()"), self, QtCore.SLOT("reject()"))

            self.exec()

        def process(self):
            plt.close('all')
            global doc
            doc = Metashape.app.document
            global app, searchRadius, maxPlaneResid, tiePoints, keyPoints, dsImageMatching, dsDenseMatching
            print("Script started...")
            self.btnP1.setDisabled(True)
            self.btnQuit.setDisabled(True)
            self.pBar.setValue(0)

            searchRadius = self.intEdt_searchRadius.text()
            if not searchRadius.isdigit():
                searchRadius = 30 # default value for the search radius
            else:
                searchRadius = int(searchRadius)

            maxPlaneResid = self.intEdt_maxPlaneResid.text()
            if not maxPlaneResid.isdigit():
                maxPlaneResid = 0.1 # default value for the search radius
            else:
                maxPlaneResid = int(maxPlaneResid)/100

            tiePoints = self.intEdt_tiePoints.text()
            if not tiePoints.isdigit():
                tiePoints = 3000 # default value for the search radius
            else:
                tiePoints = int(tiePoints)

            keyPoints = self.intEdt_keyPoints.text()
            if not keyPoints.isdigit():
                keyPoints = 25000 # default value for the search radius
            else:
                keyPoints = int(keyPoints)

            print(searchRadius)
            print(maxPlaneResid)
            print(keyPoints)
            print(tiePoints)

            selected = False
            if self.radioBtn_high.isChecked():
                dsImageMatching = 1
                dsDenseMatching = 2
            elif self.radioBtn_medium.isChecked():
                dsImageMatching = 2
                dsDenseMatching = 4
            elif self.radioBtn_slow.isChecked():
                dsImageMatching = 4
                dsDenseMatching = 16


            self.btnP1.setDisabled(False)
            self.btnQuit.setDisabled(False)
            runAll()
            return 1


    def custom_process():
        global app
        app = QtWidgets.QApplication.instance()
        parent = app.activeWindow()
        dlg = ProgDlg(parent)




    label = "AUTOGRAF"
    metashapeApp.removeMenuItem("AUTOGRAF")
    metashapeApp.addMenuItem("AUTOGRAF/1. Choose graffito directory to be processed", getPath)
    metashapeApp.addMenuItem('AUTOGRAF/-----------------------------------', dummyFunc)
    metashapeApp.addMenuItem("AUTOGRAF/2. Run", custom_process)
    metashapeApp.addMenuItem('AUTOGRAF/-----------------------------------', dummyFunc)
    metashapeApp.addMenuItem("AUTOGRAF/About AUTOGRAF", about)

    print("To execute this script press {}".format(label))






    #app.removeMenuItem("AUTOGRAF")
    #app.addMenuItem("AUTOGRAF/1. Choose graffito directory to be processed", getPath)
    #app.addMenuItem('AUTOGRAF/-----------------------------------', dummyFunc)
    #app.addMenuItem("AUTOGRAF/2. Run", runAll)
    #app.addMenuItem('AUTOGRAF/-----------------------------------', dummyFunc)
    #app.addMenuItem("AUTOGRAF/About AUTOGRAF", about)
    #
    #print('---------------------')
    #print("AUTOGRAF was added to Metashape")
    #print('---------------------')

