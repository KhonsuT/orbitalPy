from sgp4.api import Satrec
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import TEME, GCRS, EarthLocation
from astropy.coordinates.erfa_astrom import erfa

### TLEPropagator - class


class TLEPropagator:

    def __init__(self, line1, line2):
        # validate
        self.line1 = line1
        self.line2 = line2
        self.sat = Satrec.twoline2rv(self.line1, self.line2)

    def propagate(self, start_time, end_time, dt):
        t_span = np.arange(start_time, end_time, dt)
        p_res, v_res = [[], [], []], [[], [], []]
        for t in t_span:
            _, p, v = self.get_position_eci(t)
            for i in range(len(p_res)):
                p_res[i].append(p[i])
            for i in range(len(v_res)):
                v_res[i].append(v[i])
        return p_res, v_res, t_span

    def get_position_eci(self, time):
        try:
            if isinstance(time, pd.Timestamp):
                julian_date = time.to_julian_date()
            julian_date = time
            jd = math.trunc(julian_date)
            fr = julian_date - jd
            error_code, position, velocity = self.sat.sgp4(jd, fr)

            if error_code != 0:
                print("Errors in Propagation")
            return error_code, position, velocity
        except Exception as e:
            print(f"Encounter an Error of type: {e}")

    def get_position_geodetic(self, time):
        try:
            _, position, _ = self.get_position_eci(time)
            ecef_x, ecef_y, ecef_z = self.eci_to_ecef(position, time)
            lat, lon, alt = self.ecef_to_geodetic(ecef_x, ecef_y, ecef_z)
            return lat, lon, alt
        except Exception as e:
            print(f"Encounter an Error of type: {e}")

    def eci_to_ecef(self, eci, dt):
        def gmst_from_datetime(dt):
            time = Time(dt, format="jd", scale="utc")
            gmst = time.sidereal_time("mean", "greenwich").radian
            return gmst

        x_eci, y_eci, z_eci = eci
        gmst = gmst_from_datetime(dt)

        cos_theta = np.cos(gmst)
        sin_theta = np.sin(gmst)

        x_ecef = cos_theta * x_eci + sin_theta * y_eci
        y_ecef = -sin_theta * x_eci + cos_theta * y_eci
        z_ecef = z_eci

        return x_ecef, y_ecef, z_ecef

    def ecef_to_geodetic(self, x, y, z):
        a = 6378.137
        f = 1 / 298.257223563
        e2 = f * (2 - f)

        lon = np.arctan2(y, x)
        r = np.sqrt(x**2 + y**2)
        lat = np.arctan2(z, r)

        lat_prev = 0
        while abs(lat - lat_prev) > 1e-10:
            lat_prev = lat
            N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
            alt = r / np.cos(lat) - N
            lat = np.arctan2(z, r * (1 - e2 * N / (N + alt)))

        N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        alt = r / np.cos(lat) - N

        return np.degrees(lat), np.degrees(lon), alt

    def all_passes(self, ground, start_time, end_time, threshold):
        # ground station accept two types of input (eci) or (ecef)
        # ground = [lat, lon, alt] or Ground object(contains keys lat, lon, alt)
        # threshold % error i.e. 0.1 = 10% error 0.01 = 1% error
        if isinstance(ground, list):
            lat_g = ground[0]
            lon_g = ground[1]
            alt_g = ground[2]
        else:
            lat_g = ground.lat
            lon_g = ground.lon
            alt_g = ground.alt

        passes = []

        for t in np.linspace(start_time, end_time, 1000):
            lat, lon, alt = self.get_position_geodetic(t)

            lat_error = abs(lat - lat_g)
            lon_error = abs(lon - lon_g)
            alt_error = abs(alt - alt_g)

            lat_margin = threshold * abs(lat_g)
            lon_margin = threshold * abs(lon_g)
            alt_margin = threshold * abs(alt_g)
            if alt_g == 0:
                if lat_error <= lat_margin and lon_error <= lon_margin:
                    passes.append(
                        (Time(t, format="jd", scale="utc").to_datetime(), lat, lon, alt)
                    )
            else:
                if (
                    lat_error <= lat_margin
                    and lon_error <= lon_margin
                    and alt_error <= alt_margin
                ):
                    passes.append(
                        (Time(t, format="jd", scale="utc").to_datetime(), lat, lon, alt)
                    )
        return passes

    def simulate(self, end_time, start_time=datetime.now(), dt=0.001):
        # simulate functions have two modes:
        # 1. from now to future days
        # 2. from start to end
        X, _, _ = self.propagate(start_time, end_time, dt)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot(X[0], X[1], X[2], "k")
        ax.set_xlabel("x(km)")
        ax.set_ylabel("y(km)")
        ax.set_zlabel("z(km)")
        ax.set_title(f"TLE Propagation from {start_time} to {end_time}")
        plt.show()
