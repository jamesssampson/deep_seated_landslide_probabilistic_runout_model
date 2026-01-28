from .probabilistic_landslide_runout_plugin import ProbabilisticLandslideRunoutPlugin

def classFactory(iface):
    return ProbabilisticLandslideRunoutPlugin(iface)
