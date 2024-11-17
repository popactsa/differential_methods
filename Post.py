import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

step = 99
file_name = "data/" + str(step) + ".csv"
data = pd.read_csv(file_name, sep = ' ', header=None)

# print(data)

plt.plot(data.iloc[:, 0], data.iloc[:, 1], 'o-', lw = 1, markersize = 5, label = "плотность")
plt.title(f'Step = {step}')
plt.legend()
plt.show()