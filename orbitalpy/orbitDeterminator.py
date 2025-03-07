import numpy as np
from .TLEPropagator import TLEPropagator
from .FakeGPS import FakeGPS
from scipy.optimize import fmin

## Orbit Determination and update using TLE and GPS Tracking
# Combining Kalman filter with previous TLE prediction with active GPS data


class OrbitDeterminator:
    def __init__(
        self, tle_line1, tle_line2, initial_state, gps=None, state_estimator=None
    ):
        self.tle = TLEPropagator(
            line1=tle_line1, line2=tle_line2
        )  # propagator used for model prediction
        self.gps = (
            gps if gps else FakeGPS(tle_line1, tle_line2)
        )  # gps used for measurement data
        self.A = np.eye(6)  # default prediction matrix
        self.P = np.eye(6)  # default error covariance matrix
        self.H = np.eye(3, 6)  # measurement Matrix
        self.Q = np.eye(6)  # process noise
        self.R = np.eye(3)  # measurement noise
        self.X = initial_state  # initial state
        self.state_estimator = (
            state_estimator if state_estimator else self.__tle_state_estimation
        )  # pass in to use custom model

    def kalman_filter_predict_update(self, state, A, P, H, Z, R, Q):
        Xp = np.dot(A, state)

        Pp = np.dot(np.dot(A, P), A.T) + Q

        K = np.dot(np.dot(Pp, H.T), np.linalg.inv(np.dot(np.dot(H, Pp), H.T) + R))

        X = Xp + np.dot(K, (Z - np.dot(H, Xp)))
        P = np.dot(np.eye(len(state)) - np.dot(K, H), Pp)

        return X, P

    def __tle_state_estimation(self, time, dt):
        # tle state estimator
        _, pos1, vel1 = self.tle.get_position_eci(time)
        _, pos2, vel2 = self.tle.get_position_eci(time + dt)
        A = np.zeros((6, 6))
        A[:3, :3] = np.eye(3)
        A[:3, 3:] = np.eye(3) * dt
        A[3:, :3] = [(pos1[i] - pos2[i]) / dt for i in range(len(pos1))]
        A[3:, 3:] = np.eye(3)

        return A

    def __run_kalman(self, time, dt):
        # kalman iteration
        _, r, v = self.gps.get_gps_data(time)
        Z = r
        self.A = self.state_estimator(time, dt)
        self.X, self.P = self.kalman_filter_predict_update(
            self.X, self.A, self.P, self.H, Z, self.R, self.Q
        )
        return self.X, self.P, Z

    def determine(self, start_time, end_time, dt):
        # Generate position and velocity data over a period of time(julian time)
        t_span = np.arange(start_time, end_time, dt)
        X, P, Z = [[], [], []], [[], [], []], [[], [], []]
        for t in t_span:
            X_e, P_e, Z_e = self.__run_kalman(t, dt)
            for i in range(len(X)):
                X[i].append(X_e[i])
            for i in range(len(P)):
                P[i].append(P_e[i])
            for i in range(len(Z)):
                Z[i].append(Z_e[i])
        return X, P, Z
