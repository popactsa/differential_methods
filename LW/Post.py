import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file_name = "data/100.csv"
data = pd.read_csv(file_name, sep = ' ', header=None)

# print(data)

plt.plot(data.iloc[:, 0], data.iloc[:, 1])
plt.show()