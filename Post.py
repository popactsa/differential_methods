import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

step = 99
file_name = "data/" + str(step) + ".csv"
data = pd.read_csv(file_name, sep = ' ', header=None)

# print(data)

fig, ax = plt.subplots(1, 2, figsize = (15, 10), dpi = 350);
ax[0].plot(data.iloc[:, 0], data.iloc[:, 1], 'o-', lw = 1, markersize = 5, label = "плотность")
ax[0].title.set_text(f'Step = {step}')
ax[0].legend()
ax[0].grid()
ax[1].plot(data.iloc[:, 0], data.iloc[:, 2], 'o-', lw = 1, markersize = 5, label = "давление")
ax[1].title.set_text(f'Step = {step}')
ax[1].legend()
ax[1].grid()
plt.savefig("graph")
#plt.show()
