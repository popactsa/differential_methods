import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

step = 49
file_name = "data/" + str(step) + ".csv"
file_sol = "data/exact_solution.csv"
data = pd.read_csv(file_name, sep = ' ', header=None)
data1 = pd.read_csv(file_sol, sep = ' ', header=None)

# print(data)

fig, ax = plt.subplots(1, 3, figsize = (9*1.5, 6), dpi = 200);

ax[0].plot(data.iloc[:, 0], data.iloc[:, 1], 'o-', lw = 1, markersize = 3, label = "плотность")
ax[0].plot(data1.iloc[:, 0], data1.iloc[:, 1], '-', c = 'k', lw = 1, markersize = 3, label = "точное решение")
ax[0].title.set_text(f'Step = {step}')
ax[0].legend(fontsize=6.0)
ax[0].grid()

ax[1].plot(data.iloc[:, 0], data.iloc[:, 2], 'o-', lw = 1, markersize = 3, label = "скорость")
ax[1].plot(data1.iloc[:, 0], data1.iloc[:, 2], '-', c = 'k', lw = 1, markersize = 3, label = "точное решение")
ax[1].title.set_text(f'Step = {step}')
ax[1].legend(fontsize=6.0)
ax[1].grid()

ax[2].plot(data.iloc[:, 0], data.iloc[:, 3], 'o-', lw = 1, markersize = 3, label = "давление")
ax[2].plot(data1.iloc[:, 0], data1.iloc[:, 3], '-', c = 'k', lw = 1, markersize = 3, label = "точное решение")
ax[2].title.set_text(f'Step = {step}')
ax[2].legend(fontsize=6.0)
ax[2].grid()

plt.savefig("graph")
# plt.show()
