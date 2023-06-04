
from physt import histogram, binnings
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(42)

heights1 = np.random.normal(169, 10, 100000)
heights2 = np.random.normal(180, 6, 100000)
numbers = np.random.rand(100000)

plt.rcParams.update({'font.size': 18})
figure, axis = plt.subplots(1, 2, figsize=(10, 4))
# bins2 = binning.quantile_bins(heights1, 40)
hist2 = histogram(heights1, "quantile", bin_count=20)
hist2.plot(density=True, ax=axis[1])
axis[1].set_title("Equal depth binning")
# axis[1].grid(which='both',axis='x', linewidth=0.25, color='Black')
# histogram(heights1, 10).plot(alpha=0.3, density=True, ax=axis[1], color="green", label="Equal spaced")
# axis[1].legend(loc=2)
plt.show()
