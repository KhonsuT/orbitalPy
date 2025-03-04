from sgp4.api import Satrec
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import math
### TLEPropagator - class

tle_line1 = "1 25544U 98067A   25062.66313727  .00016275  00000-0  29551-3 0  9996"
tle_line2 = "2 25544  51.6371 109.6890 0005902 339.6584 166.5761 15.49702064498799"

class TLEPropagator:

    def __init__(self, line1, line2):
        # validate
        self.line1 = line1
        self.line2 = line2
        self.sat = Satrec.twoline2rv(self.line1, self.line2)

    def get_position_eci(self, time):
        julian_date = time.to_julian_date()
        jd = math.trunc(julian_date)
        fr = julian_date - jd
        error_code, position, velocity = self.sat.sgp4(jd,fr)

        if error_code != 0:
            print("Errors in Propagation")
        return error_code, position, velocity
    
    def get_position_ecef(self, time):
        error_code, position, velocity = self.get_position_eci(time)
        lat, lon, alt = self.eci_to_ecef(position)
        return lat, lon, alt

    def eci_to_ecef(self,*args):
        if len(args) == 1:
            x = args[0][0]
            y = args[0][1]
            z = args[0][2]
        elif len(args) == 3:
            x = args[0]
            y = args[1]
            z = args[2]
        else:
            raise ValueError("Invalid input for eci_to_ecef conversion, must be one of [position(x,y,z),[x,y,z]]")
        a = 6378.137
        deg2rad = np.pi / 180
        r = np.sqrt(x**2 + y**2)
        lon = np.arctan2(y,x) / deg2rad
        lat = np.arctan2(z,r) / deg2rad
        alt = np.sqrt(x**2 + y**2 + z**2) - a

        return lat, lon, alt

    def all_passes(self, ground, start_time, end_time):
        ### Function that returns a list of passes based on inputed time and ground ref
        ## accuracy to 5 decimal points
        

    
sat = TLEPropagator(tle_line1,tle_line2)

print(sat.get_position_ecef(pd.to_datetime('2023-08-15 12:30:00')))