## Authors: Colin Mazengarb(Original Determanistic Model plugin is based on)
##          James Sampson(Integration of original model into a probabalistic one and formation of plugin)
## Purpose of Script: A script that models the runout of deep seated landslides
## The script uses an alpha value derived from the geology vector layer to control.
# Runout is controlled by three limiters:
# 1) the Conefall method where runout is not so strongly controlled by local topography.
# 2) total travel distance:source distance ratio of 3.5 to prevent unreasonably long runouts (developed for the NW Coast of Tasmania),  particularly where thin basaltic units (potentially susceptible to failure) overlie non-susceptible basement units that are deeply engorged.
# 3) An arbitrary length limiter as a further control on the runout to prevent unreasonable runouts
## All input raster file formats are to be tifs, all inputs are aligned to the DEM within the plugin
## "source" = This is created through a query of a slope grid; shall consist of sources areas = 1, non source areas = 0 or null. For testing this is taken from MRT
## "dem" = the name of the digital elevation model; any floating point values or null.
## "Geology" = Alpha angles derived from tasmanian geology layers and values from published MRT litrature
#Import nessasary libaries
from qgis.core import (
    Qgis,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterNumber,
    QgsRasterLayer,
    QgsProcessingLayerPostProcessorInterface,
    QgsProcessingContext,
    QgsProcessingParameterBoolean
)

from qgis.PyQt.QtCore import QCoreApplication
import processing
import math
import time
import numpy as np
import os


class ProbRunoutPostProcessor(QgsProcessingLayerPostProcessorInterface):
    def __init__(self, style_path):
        super().__init__()
        self.style_path = style_path

    def postProcessLayer(self, layer, context, feedback):
        if isinstance(layer, QgsRasterLayer) and os.path.exists(self.style_path):
            layer.loadNamedStyle(self.style_path)
            layer.triggerRepaint()


class LandslideRunoutAlgorithmMC(QgsProcessingAlgorithm):

    SRC = 'SRC'
    DEM = 'DEM'
    FRIC = 'FRIC'
    OUTPUT = 'OUTPUT'

    FAN_ANGLE = 'FAN_ANGLE'
    FAN_INCREMENT = 'FAN_INCREMENT'
    RUNOUT_LIMITER = 'RUNOUT_LIMITER'
    RATIO_LIMITER = 'RATIO_LIMITER'

    FRIC_COEFF = 'FRIC_COEFF'
    SLOPE_COEFF = 'SLOPE_COEFF'
    ALPHA_MIN = 'ALPHA_MIN'
    ALPHA_MAX = 'ALPHA_MAX'

    N_REALIZ = 'N_REALIZ'  # Allows the user to change the number of Monte Carlo runs
    MODE = 'MODE'          # Allows user to choose probabilistic vs deterministic

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def name(self):
        return 'landslide_runout_alpha_mc_fric'

    def displayName(self):
        return self.tr('Landslide runout (for deep seated Tasmanian landslides)')

    def group(self):
        return self.tr('Landslides')

    def groupId(self):
        return 'landslides'

    def createInstance(self):
        return LandslideRunoutAlgorithmMC()

    def initAlgorithm(self, config=None):

        self.addParameter(QgsProcessingParameterRasterLayer(
            self.SRC, self.tr('Source raster (0/1)')))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.DEM, self.tr('DEM')))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.FRIC, self.tr('Friction raster')))
        self.addParameter(QgsProcessingParameterRasterDestination(
            self.OUTPUT, self.tr('Runout impact probability (0–1)')))

        self.addParameter(QgsProcessingParameterNumber(
            self.FAN_ANGLE, self.tr('Fan angle (deg)'),
            QgsProcessingParameterNumber.Double, 20))
        self.addParameter(QgsProcessingParameterNumber(
            self.FAN_INCREMENT, self.tr('Fan increment (deg)'),
            QgsProcessingParameterNumber.Double, 10))
        self.addParameter(QgsProcessingParameterNumber(
            self.RUNOUT_LIMITER, self.tr('Runout limit (m)'),
            QgsProcessingParameterNumber.Double, 50))
        self.addParameter(QgsProcessingParameterNumber(
            self.RATIO_LIMITER, self.tr('Runout:source ratio limit'),
            QgsProcessingParameterNumber.Double, 3.5))

        self.addParameter(QgsProcessingParameterNumber(
            self.FRIC_COEFF, self.tr('Friction coefficient'),
            QgsProcessingParameterNumber.Double, 1.0))
        self.addParameter(QgsProcessingParameterNumber(
            self.SLOPE_COEFF, self.tr('Slope coefficient'),
            QgsProcessingParameterNumber.Double, 0.5))
        self.addParameter(QgsProcessingParameterNumber(
            self.ALPHA_MIN, self.tr('Min alpha (deg)'),
            QgsProcessingParameterNumber.Double, 10.0))
        self.addParameter(QgsProcessingParameterNumber(
            self.ALPHA_MAX, self.tr('Max alpha (deg)'),
            QgsProcessingParameterNumber.Double, 45.0))

        # Monte Carlo iterations as integer parameter
        self.addParameter(QgsProcessingParameterNumber(
            self.N_REALIZ,
            self.tr('Monte Carlo iterations'),
            Qgis.ProcessingNumberParameterType.Integer,
            100,    # default
            False,  # not optional
            1,      # min
            10000   # max
        ))

        # Allow choice between probabilistic or deterministic
        self.addParameter(QgsProcessingParameterBoolean(
            self.MODE,
            self.tr('Use Probabilistic Model (leave unchecked to use Determanistic Model)'),
            defaultValue=True
        ))

    def processAlgorithm(self, parameters, context, feedback):

        from osgeo import gdal
        from osgeo.gdalconst import GDT_Float32

        gdal.SetConfigOption('CHECK_DISK_FREE_SPACE', 'FALSE')

        start_time = time.time()

        #Checks the rasters taken from the
        src_layer = self.parameterAsRasterLayer(parameters, self.SRC, context)
        dem_layer = self.parameterAsRasterLayer(parameters, self.DEM, context)
        fric_layer = self.parameterAsRasterLayer(parameters, self.FRIC, context)

        for lyr, name in [(src_layer, 'SRC'), (dem_layer, 'DEM'), (fric_layer, 'FRIC')]:
            if lyr is None or not lyr.isValid():
                raise QgsProcessingException(f'{name} raster is invalid')

        src_path = src_layer.dataProvider().dataSourceUri()
        dem_path = dem_layer.dataProvider().dataSourceUri()
        fric_path = fric_layer.dataProvider().dataSourceUri()

        out_path = self.parameterAsOutputLayer(parameters, self.OUTPUT, context)

        #Numeric parameters to feed into script
        FanAngle = self.parameterAsDouble(parameters, self.FAN_ANGLE, context)
        FanIncrement = self.parameterAsDouble(parameters, self.FAN_INCREMENT, context)
        RunoutLimiter = self.parameterAsDouble(parameters, self.RUNOUT_LIMITER, context)
        RunoutRatioLimiter = self.parameterAsDouble(parameters, self.RATIO_LIMITER, context)

        fric_coeff = self.parameterAsDouble(parameters, self.FRIC_COEFF, context)
        slope_coeff = self.parameterAsDouble(parameters, self.SLOPE_COEFF, context)
        alpha_min = self.parameterAsDouble(parameters, self.ALPHA_MIN, context)
        alpha_max = self.parameterAsDouble(parameters, self.ALPHA_MAX, context)

        # Monte Carlo iterations from parameter
        N = self.parameterAsInt(parameters, self.N_REALIZ, context)

        # Mode: probabilistic vs deterministic
        use_mc = self.parameterAsBool(parameters, self.MODE, context)

        #Takes DEM and creates slope and aspect and puts them into temp output for script use
        slope_path = processing.run(
            "qgis:slope",
            {'INPUT': dem_layer, 'Z_FACTOR': 1.0, 'OUTPUT': 'TEMPORARY_OUTPUT'},
            context=context, feedback=feedback, is_child_algorithm=True
        )['OUTPUT']

        aspect_path = processing.run(
            "qgis:aspect",
            {'INPUT': dem_layer, 'Z_FACTOR': 1.0, 'OUTPUT': 'TEMPORARY_OUTPUT'},
            context=context, feedback=feedback, is_child_algorithm=True
        )['OUTPUT']

        #Opens gdal(Geospatial Data Abstraction Libary)
        srcDS = gdal.Open(src_path)
        demDS = gdal.Open(dem_path)
        fricDS = gdal.Open(fric_path)
        slopeDS = gdal.Open(slope_path)
        aspectDS = gdal.Open(aspect_path)

        if any(ds is None for ds in [srcDS, demDS, fricDS, slopeDS, aspectDS]):
            raise QgsProcessingException('Failed to open one or more rasters')

        transform = demDS.GetGeoTransform()
        pixel = transform[1]

        rows = demDS.RasterYSize
        cols = demDS.RasterXSize

        src = srcDS.GetRasterBand(1).ReadAsArray().astype(np.uint8)
        dem = demDS.GetRasterBand(1).ReadAsArray().astype(np.float32)
        fric_base = fricDS.GetRasterBand(1).ReadAsArray().astype(np.float32)
        slope = slopeDS.GetRasterBand(1).ReadAsArray().astype(np.float32)
        aspect = aspectDS.GetRasterBand(1).ReadAsArray().astype(np.float32)

        srcDS = None
        demDS = None
        fricDS = None
        slopeDS = None
        aspectDS = None

        total_sources = int((src == 1).sum())
        feedback.pushInfo(f'Source cells: {total_sources}')

        #Settings for Monte Carlo. This can be changed to increase the distance from 'truth' ie can increase the range of values
        fric_sigma_factor = 0.1
        fric_min = np.nanpercentile(fric_base, 1) - 3.0
        fric_max = np.nanpercentile(fric_base, 99) + 3.0

        impact_counts = np.zeros((rows, cols), np.int32)

        if use_mc:
            #Probabilistic mode
            for n in range(N):
                if feedback.isCanceled():
                    break

                sigma = fric_sigma_factor * np.maximum(np.abs(fric_base), 1e-3)
                noise = np.random.normal(
                    loc=0.0,
                    scale=sigma,
                    size=fric_base.shape
                ).astype(np.float32)

                fric = fric_base + noise
                fric = np.clip(fric, fric_min, fric_max)

                alpha = np.clip(
                    fric_coeff * fric - slope_coeff * slope,
                    alpha_min, alpha_max
                )

                runout_n = self._deterministic_runout(
                    src, dem, slope, aspect, alpha,
                    FanAngle, FanIncrement,
                    RunoutLimiter, RunoutRatioLimiter,
                    pixel, rows, cols,
                    feedback, N, n
                )

                impact_counts += (runout_n > 0).astype(np.int32)
                feedback.setProgress(int((n + 1) / N * 100))
        else:
            #Determanistic mode, no use of MC and no noise added
            fric = fric_base.copy()

            alpha = np.clip(
                fric_coeff * fric - slope_coeff * slope,
                alpha_min, alpha_max
            )

            runout_n = self._deterministic_runout(
                src, dem, slope, aspect, alpha,
                FanAngle, FanIncrement,
                RunoutLimiter, RunoutRatioLimiter,
                pixel, rows, cols,
                feedback, 1, 0
            )

            impact_counts = (runout_n > 0).astype(np.int32)
            N = 1  # so impact_prob is 0/1 in deterministic mode

        #Probability
        impact_prob = impact_counts.astype(np.float32) / float(N)

        #Writes output
        outDS = gdal.GetDriverByName('GTiff').Create(
            out_path, cols, rows, 1, GDT_Float32,
            options=['COMPRESS=LZW', 'TILED=YES']
        )

        band = outDS.GetRasterBand(1)
        band.SetNoDataValue(-9999)
        band.WriteArray(impact_prob)

        outDS.SetGeoTransform(transform)
        outDS.SetProjection(src_layer.crs().toWkt())
        outDS = None

        feedback.pushInfo(
            f'Completed {N} simulations in {round((time.time() - start_time) / 60, 2)} minutes'
        )

        return {self.OUTPUT: out_path}

    def _deterministic_runout(self, src, dem, slope, aspect, alpha,
                              FanAngle, FanIncrement,
                              RunoutLimiter, RunoutRatioLimiter,
                              pixel, rows, cols,
                              feedback, N, n):
        #Core Deterministic runout used for both models
        runout_n = np.zeros((rows, cols), np.float32)

        for i in range(rows):
            if feedback.isCanceled():
                break

            for j in range(cols):
                if src[i, j] != 1:
                    continue

                runout_n[i, j] = 1.0

                elev_src = dem[i, j]
                alpha_angle = alpha[i, j]

                for fan in range(-int(FanAngle), int(FanAngle), int(FanIncrement)):
                    x, y = j, i
                    L = 0.0
                    src_len = 0.0

                    direction = (aspect[i, j] + fan) % 360
                    dx = math.sin(math.radians(direction))
                    dy = -math.cos(math.radians(direction))

                    for _ in range(1000):
                        x += dx
                        y += dy
                        L += pixel

                        rx, ry = int(round(x)), int(round(y))
                        if rx < 0 or ry < 0 or rx >= cols or ry >= rows:
                            break

                        A = math.degrees(
                            math.atan((elev_src - dem[ry, rx]) / L)
                        )

                        if src[ry, rx] != 1:
                            src_len += pixel

                        if (
                            A < alpha_angle or
                            L > RunoutLimiter or
                            (src_len > 0 and L / src_len > RunoutRatioLimiter)
                        ):
                            break

                        runout_n[ry, rx] = max(runout_n[ry, rx], L)
            if N > 0:
                progress = int(((n + (i + 1) / rows) / N) * 100)
                feedback.setProgress(progress)

        return runout_n

    def postProcessAlgorithm(self, context: QgsProcessingContext, feedback):
        style_path = os.path.join(
            os.path.dirname(__file__), '..', 'styles', 'probability_25_50_90.qml'
        )

        if not os.path.exists(style_path):
            return {}

        layers = context.layersToLoadOnCompletion()
        # self.OUTPUT is the key that returned in processAlgorithm
        if self.OUTPUT in layers:
            details = layers[self.OUTPUT]
            pp = ProbRunoutPostProcessor(style_path)
            details.setPostProcessor(pp)

        return {}
