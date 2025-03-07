import numpy as np
from .TLEPropagator import TLEPropagator


class FakeGPS:
    def __init__(self, tle_line1, tle_line2):
        # Fake gps to simulate noisy gps input data
        self.GPSPropagator = TLEPropagator(line1=tle_line1, line2=tle_line2)
        self.mean = 0
        self.std_dev = 25

    def get_gps_data(self, time):
        # input julian time
        e, r, v = self.GPSPropagator.get_position_eci(time)
        r = [val + np.random.normal(self.mean, self.std_dev) for val in r]
        v = [val + np.random.normal(self.mean, self.std_dev) for val in v]
        return e, r, v
