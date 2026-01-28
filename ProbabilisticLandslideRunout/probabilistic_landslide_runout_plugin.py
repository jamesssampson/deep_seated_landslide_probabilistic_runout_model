from qgis.core import QgsApplication
from .processing_provider import ProbabilisticLandslideRunoutProvider

class ProbabilisticLandslideRunoutPlugin:

    def __init__(self, iface):
        self.iface = iface
        self.provider = None

    def initGui(self):
        self.provider = ProbabilisticLandslideRunoutProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def unload(self):
        if self.provider is not None:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None
