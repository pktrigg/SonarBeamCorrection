import SonarBeamCorrection 
from scipy import interpolate
# import time
import numpy as np
import matplotlib.pyplot as plt
import math

def main():
    """here we round trip the ping to beam correction algorithm to prove we can flip between a ping, angular response curve and back again """

    pingIn = np.arange(1,4096, dtype=np.int  )
    slantRange =1000
    altitude = 5
    ARCStartAngle = 0.1
    ARCEndAngle = 90
    ARCAngleInterval = 0.1

    ARC = SonarBeamCorrection.pingToAngularResponse(pingIn, slantRange, altitude, ARCStartAngle, ARCEndAngle, ARCAngleInterval)
    #now convert back to a ping
    pingOut = SonarBeamCorrection.angularResponseToPing(ARC, slantRange, altitude, pingIn.size, ARCStartAngle, ARCEndAngle, ARCAngleInterval)

    plt.xlabel('sample')
    plt.ylabel('intensity')
    plt.grid(True)
    plt.plot(pingIn, label="PingIn", color='b', linewidth=0.1)
    # plt.pause(1)
    # plt.plot (ARC, marker='x', label="ARC", color='r', linewidth=0.5)
    

    plt.plot(pingOut, marker='x', label="Ping Out", color='green', linewidth=0.1)
    # for i in range(0,len(pingOut)):
    #     if i < 50:
    #         plt.annotate('(%d)' % i, (i, pingOut[i]))
    plt.legend()
    plt.show()





if __name__ == "__main__":
    main()

