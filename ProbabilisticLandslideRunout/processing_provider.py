from qgis.core import QgsProcessingProvider
from .algorithms.landslide_runout_alpha_mc_fric import LandslideRunoutAlgorithmMC

class ProbabilisticLandslideRunoutProvider(QgsProcessingProvider):

    def loadAlgorithms(self):
        self.addAlgorithm(LandslideRunoutAlgorithmMC())

    def id(self):
        return 'probabilisticlandsliderunout'

    def name(self):
        return self.tr('Landslides')

    def longName(self):
        return self.name()
