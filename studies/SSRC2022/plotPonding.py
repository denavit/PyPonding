import matplotlib.pyplot as plt
import numpy as np
import sys


def plotPonding(heights):

  legends = []
  for height in heights:
    data = np.loadtxt(f"disp-{height}in.txt")
    plt.plot(data[:, 0], data[:, 1])
    legends.append(f"{height} in")

  plt.xlabel('Position along the length of beam (ft)')
  plt.ylabel('Displacement (in)')
  plt.legend(legends)
  plt.show()


if __name__ == "__main__":
  if len(sys.argv) < 2:
    print("Required:", sys.argv[0], "waterHeightInch")
    exit()

  waterHeightInch = []
  try:
    for i in range(1, len(sys.argv)):
      waterHeightInch.append(float(sys.argv[i]))
  except:
    print(f"Invalid waterHeightInch")

  plotPonding(waterHeightInch)
