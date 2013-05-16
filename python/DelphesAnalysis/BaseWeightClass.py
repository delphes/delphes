# Every weight class must contain a weight method, taking as first argument the event.
# Additional methods and arguments are allowed.
# Note that weight classes do not have to derive from BaseWeightClass.

class BaseWeightClass:
  """A class to reweight AnalysisEvents"""
  def weight( self, event):
     """Lepton eff weight"""
     raise NotImplementedError

