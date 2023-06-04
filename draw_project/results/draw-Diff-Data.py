import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker

if __name__ == '__main__':
    SimQuantileNum = 99
    SeriesNum,SeriesPerDataset = 45,15

    Diff=np.zeros((SeriesNum, SeriesNum))
    print(Diff)

    with open('Data-Diff.txt', encoding='utf-8') as file_obj:
        contents = file_obj.read()
        # lines=contents.splitlines()
        mId=0
        for line in contents.splitlines():
            values = line.split('\t')
            # print(values[0])
            nId=0
            for value in values:
                value=float(value)
                Diff[mId, nId]=value
                nId+=1
            mId+=1
    print(Diff)
    plt.close()
    plt.rcParams.update({'font.size': 18})
    ax=plt.subplot(111)
    X,Y=np.meshgrid(np.arange(SeriesNum),np.arange(SeriesNum))
    print(X)
    print(Y)
    pcm=ax.pcolor(Diff, norm=colors.PowerNorm(gamma=0.3), cmap='Blues', shading='auto')
    ax.grid(which='both',linewidth=0.25,color='Black')
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1.0))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(1.0))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(SeriesPerDataset))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(SeriesPerDataset))
    cb=plt.colorbar(pcm)
    cb.set_ticks([0,0.5])
    cb.update_ticks()
    plt.show()
