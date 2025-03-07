from orbitalpy.orbitDeterminator import OrbitDeterminator
from orbitalpy.TLEPropagator import TLEPropagator
from astropy.time import Time
from datetime import datetime, timedelta
from datetime import timezone
import numpy as np
import matplotlib.pyplot as plt
from orbitalpy.FakeGPS import FakeGPS


def main():
    #ISS tle data from 08-22-24
    tle_line1 = "1 25544U 98067A   24235.16000531  .00028314  00000+0  49252-3 0  9994"
    tle_line2 = "2 25544  51.6402 348.4828 0004893 236.4095 123.6428 15.50549723468780"
    tle = TLEPropagator(tle_line1, tle_line2)
    now = datetime.now(timezone.utc)
    julian_time = Time(now).jd
    _, position, velocity = tle.get_position_eci(julian_time)
    obd = OrbitDeterminator(tle_line1, tle_line2, np.hstack([position, velocity]))
    #ISS tle data from 08-22-24
    obd.gps = FakeGPS(
        "1 25544U 98067A   25065.79262731  .00012468  00000-0  22778-3 0  9995",
        "2 25544  51.6359  94.1869 0006272 350.9448 346.3820 15.49773953499278",
    )
    future_date = now + timedelta(days=1)
    future_julian_time = Time(future_date).jd
    true_state = [[], [], [], [], [], []]
    #ISS tle data from 08-22-24
    tle_true = TLEPropagator(
        "1 25544U 98067A   25065.79262731  .00012468  00000-0  22778-3 0  9995",
        "2 25544  51.6359  94.1869 0006272 350.9448 346.3820 15.49773953499278",
    )

    dt = 0.01
    t_span = np.arange(julian_time, future_julian_time, dt)
    X, _, _ = obd.determine(julian_time, future_julian_time, dt)
    for t in t_span:
        _, p, v = tle_true.get_position_eci(t)
        for i in range(len(p)):
            true_state[i].append(p[i])
        for i in range(len(v)):
            true_state[i + 3].append(v[i])
    o_pos, _, _ = obd.tle.propagate(julian_time, future_julian_time, dt)
    print(f"X error: {np.sqrt(np.mean((np.array(true_state[0]) - np.array(X[0]))**2))}")
    print(f"Y error: {np.sqrt(np.mean((np.array(true_state[1]) - np.array(X[1]))**2))}")
    print(f"Z error: {np.sqrt(np.mean((np.array(true_state[2]) - np.array(X[2]))**2))}")

    print("----------------------------------------")

    print(
        f"X error: {np.sqrt(np.mean((np.array(true_state[0]) - np.array(o_pos[0]))**2))}"
    )
    print(
        f"Y error: {np.sqrt(np.mean((np.array(true_state[1]) - np.array(o_pos[1]))**2))}"
    )
    print(
        f"Z error: {np.sqrt(np.mean((np.array(true_state[2]) - np.array(o_pos[2]))**2))}"
    )
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.plot(t_span, X[0], "b")
    ax1.plot(t_span, o_pos[0], "--k")
    ax1.plot(t_span, true_state[0], "r")
    ax2.plot(t_span, X[1], "b")
    ax2.plot(t_span, o_pos[1], "--k")
    ax2.plot(t_span, true_state[1], "r")
    ax3.plot(t_span, X[2], "b")
    ax3.plot(t_span, o_pos[2], "--k")
    ax3.plot(t_span, true_state[2], "r")
    plt.legend(["Estimated State", "Original State", "True State"])
    plt.show()

if __name__ == "__main__":
    main()
