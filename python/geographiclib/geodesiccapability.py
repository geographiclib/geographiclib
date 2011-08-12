class GeodesicCapability(object):
  """Capability constants for Geodesic"""

  CAP_NONE = 0
  CAP_C1   = 1<<0
  CAP_C1p  = 1<<1
  CAP_C2   = 1<<2
  CAP_C3   = 1<<3
  CAP_C4   = 1<<4
  CAP_ALL  = 0x1F
  OUT_ALL  = 0x7F80
  NONE          = 0
  LATITUDE      = 1<<7  | CAP_NONE
  LONGITUDE     = 1<<8  | CAP_C3
  AZIMUTH       = 1<<9  | CAP_NONE
  DISTANCE      = 1<<10 | CAP_C1
  DISTANCE_IN   = 1<<11 | CAP_C1 | CAP_C1p
  REDUCEDLENGTH = 1<<12 | CAP_C1 | CAP_C2
  GEODESICSCALE = 1<<13 | CAP_C1 | CAP_C2
  AREA          = 1<<14 | CAP_C4
  ALL           = OUT_ALL| CAP_ALL

