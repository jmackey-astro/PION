
import numpy as np

# 256x128 results from nova
grid = 32768.0
time=[8.2724168e+02,2.8836392e+03,6.3653492e+03,1.2861195e+04,2.2944817e+04]
zones=[1,2,3,4,5]
steps=[7090,14488,29124,58312,118976]
update_rate = [1., 1.+1./2, 1.+1./2+1./4, 1.+1./2+1./4+1./8, 1.+1./2+1./4+1./8+1./16]

time = np.array(time)
zones = np.array(zones)
steps = np.array(steps)
udpate_rate = grid*np.array(update_rate)

zone_updates = steps*update_rate*grid
work_ratio = zone_updates/zone_updates[0]
zps = zone_updates/time
time_ratio = time/time[0]

print zones
print time
print steps
print grid*zones
print work_ratio
print time_ratio
print zps

