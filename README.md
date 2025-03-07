# This library is created integrating the SGP4 TLE Propagator Library
The primary goal is to simplify TLE propation usage in python with few added capabilities

## TLE Propagator
Note the entire library works primarily in Julian time
- A wrapper class around SGP4 with functions to get position and velocity data at T
- Generate position and velocity data in ECI and geodetic frames
- simulation using matplotlib for visualization
- propagtion method for simulating extended periods of time

## OrbitDeterminator
A class dedicated to provide estimation update on satellite position and velocity using kalman filter
- Insert GPS class
- Insert Dynamic Models
  
Example Usage: 
 - A real GPS reading(insert into OrbitDeterminator class)
 - Previous Day TLE data
 - Run in continuous interval(i.e. 1 min update) to get position and velocity estimate
 - Using estimation to update TLE accordingly(i.e. least-square fitting)

