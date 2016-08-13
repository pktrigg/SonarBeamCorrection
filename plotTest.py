import numpy as np
import matplotlib.pyplot as plt

# plt.axis([0, 1, 0, 1])
plt.ion()
plt.title("interactive test")
plt.xlabel("index")

for i in range(1000):
    y = np.random.random()
    x = np.random.random()
    scat = plt.scatter(x, y)
    # plt.draw()
    plt.pause(0.05)
    print (i)
    # scat.remove()

# while True:
#     plt.pause(0.01)