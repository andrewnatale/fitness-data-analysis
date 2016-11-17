import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
import numpy as np
from numpy import ma

class MidpointNormalize(Normalize):
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint
            resdat[resdat>0] /= abs(vmax - midpoint)
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)

        if is_scalar:
            result = result[0]
        return result

def plot_map(data, row_labels, col_labels, mpoint, show_cbar=False):
    xmin, xmax = 0, data.shape[1]
    ymin, ymax = 0, data.shape[0]

    masked_array = np.ma.array(data, mask=np.isnan(data))
    cmap = matplotlib.cm.coolwarm_r
    cmap.set_bad('black')

    max_score = max(data.flatten())
    min_score = min(data.flatten())

    fig = plt.figure()
    ax = plt.imshow(masked_array, extent=[xmin, xmax, ymin, ymax], cmap=cmap,
    norm=MidpointNormalize(midpoint=mpoint, vmin = min_score, vmax = max_score), interpolation='nearest', aspect='equal')

    # Major ticks
    plt.yticks(np.arange(0.5, data.shape[0]+0.5, 1), row_labels, fontsize=10)
    plt.xticks(np.arange(0.5, data.shape[1]+0.5, 1), col_labels, fontsize=10, rotation='vertical')
    # Minor ticks
    plt.axes().set_yticks(np.arange(ymin, ymax, 1), minor=True)
    plt.axes().set_xticks(np.arange(xmin, xmax, 1), minor=True)

    if show_cbar:
        plt.colorbar()

    plt.axes().spines['bottom'].set_color('w')
    plt.axes().spines['top'].set_color('w')
    plt.axes().spines['right'].set_color('w')
    plt.axes().spines['left'].set_color('w')
    plt.grid(which='minor', color='w', linestyle='-', linewidth=0.3)
    plt.grid(which='major', color='w', linestyle='-', linewidth=0)

    plt.tick_params(axis=u'both', which=u'both',length=0)
    plt.tight_layout()
    plt.show()

# test_array = np.random.rand(5,7)
# print test_array.shape
# col_labels = range(0,test_array.shape[1])
# row_labels = [i+0.2 for i in range(0,test_array.shape[0])]
# print col_labels, row_labels
#
# plot_map(test_array, row_labels, col_labels)
