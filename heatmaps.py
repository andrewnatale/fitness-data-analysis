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

# plot one heatmap in figure
def single_map(data, row_labels, col_labels, desc=None, use_mpoint=False, mpoint=0.0, set_limits=None, show_cbar=False):
    xmin, xmax = 0, data.shape[1]
    ymin, ymax = 0, data.shape[0]
    masked_array = np.ma.array(data, mask=np.isnan(data))

    if use_mpoint:
        cmap = matplotlib.cm.coolwarm_r
        print('using cmap: coolwarm')
    else:
        cmap = matplotlib.cm.plasma
        print('using cmap: plasma')

    cmap.set_bad('black')

    if set_limits:
        max_score, min_score = set_limits
        print('user limits specified: ', max_score, min_score)
    else:
        max_score = max(data.flatten())
        min_score = min(data.flatten())
        print('automatic limits: ', max_score, min_score)

    fig = plt.figure()
    if use_mpoint:
        print('using defined midpoint to normalize colormap')
        ax = plt.imshow(masked_array, extent=[xmin, xmax, ymin, ymax], cmap=cmap,
        norm=MidpointNormalize(midpoint=mpoint, vmin = min_score, vmax = max_score), interpolation='nearest', aspect='equal')
    else:
        print('no midpoint given, ploting unidirectional colormap')
        ax = plt.imshow(masked_array, extent=[xmin, xmax, ymin, ymax], cmap=cmap,
        vmin = min_score, vmax = max_score, interpolation='nearest', aspect='equal')

    if desc:
        plt.axes().set_title(desc)

    # Major ticks
    # b/c matplotlib is dumb, flip the row labels so they match the heatmap orientation
    row_labels.reverse()
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

# plot multiple heatmaps side by side - they will all use the same options
def new_multi_map(data_list, row_labels, col_labels, desc_list=None, use_mpoint=False, mpoint=0.0, set_limits=None, show_cbar=False):
    count_plots = len(data_list)
    if use_mpoint:
        cmap = matplotlib.cm.coolwarm
        print('using cmap: coolwarm')
    else:
        cmap = matplotlib.cm.plasma
        print('using cmap: plasma')
    cmap.set_bad('black')

    plot_id = 1
    fig = plt.figure()
    for data in data_list:
        xmin, xmax = 0, data.shape[1]
        ymin, ymax = 0, data.shape[0]

        masked_array = np.ma.array(data, mask=np.isnan(data))

        if set_limits:
            max_score, min_score = set_limits
            print('user limits specified: ', max_score, min_score)
        else:
            max_score = max(data.flatten())
            min_score = min(data.flatten())
            print('automatic limits: ', max_score, min_score)

        ax = fig.add_subplot(1,count_plots,plot_id)
        if use_mpoint:
            print('using defined midpoint to normalize colormap')
            plt.imshow(masked_array, extent=[xmin, xmax, ymin, ymax], cmap=cmap,
            norm=MidpointNormalize(midpoint=mpoint, vmin = min_score, vmax = max_score), interpolation='nearest', aspect='equal')
        else:
            print('no midpoint given, ploting unidirectional colormap')
            plt.imshow(masked_array, extent=[xmin, xmax, ymin, ymax], cmap=cmap,
            vmin = min_score, vmax = max_score, interpolation='nearest', aspect='equal')

        if desc_list:
            ax.set_title(desc_list[plot_id-1])

        # Major ticks
        # b/c matplotlib is dumb, flip the row labels so they match the heatmap orientation
        row_labels.reverse()
        plt.yticks(np.arange(0.5, data.shape[0]+0.5, 1), row_labels, fontsize=10)
        plt.xticks(np.arange(0.5, data.shape[1]+0.5, 1), col_labels, fontsize=10, rotation='vertical')
        # Minor ticks
        ax.set_yticks(np.arange(ymin, ymax, 1), minor=True)
        ax.set_xticks(np.arange(xmin, xmax, 1), minor=True)

        ax.spines['bottom'].set_color('w')
        ax.spines['top'].set_color('w')
        ax.spines['right'].set_color('w')
        ax.spines['left'].set_color('w')

        ax.grid(which='minor', color='w', linestyle='-', linewidth=0.3)
        ax.grid(which='major', color='w', linestyle='-', linewidth=0)
        plt.tick_params(axis=u'both', which=u'both',length=0)
        plot_id += 1
    if show_cbar:
        cbaxes = fig.add_axes([0.1, 0.1, 0.8, 0.03])
        cb = plt.colorbar(orientation='horizontal', cax = cbaxes)
        # plt.subplot(1,count_plots,plot_id)
        #plt.colorbar(shrink=0.7, pad=0.5)
    plt.tight_layout()
    plt.show()
