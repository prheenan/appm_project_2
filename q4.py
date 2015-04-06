# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
# need to add the utilities class. Want 'home' to be platform independent
from os.path import expanduser
home = expanduser("~")
# get the utilties directory (assume it lives in ~/utilities/python)
# but simple to change
path= home +"/utilities/python"
import sys
sys.path.append(path)
# import the patrick-specific utilities
import GenUtilities  as pGenUtil
import PlotUtilities as pPlotUtil
import CheckpointUtilities as pCheckUtil

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


def plotKMesh(length,kVals,q,ax):
    kMesh,qMesh = np.meshgrid(kVals,q)
    func = (length-kMesh+1) * qMesh**(kMesh)
    surf = ax.plot_trisurf(kMesh.flatten(), qMesh.flatten(), func.flatten(),
                           cmap=cm.hot,linewidth=0.0,antialiased=True)
    ax.set_xlabel('k-mer size')
    ax.set_ylabel('q, maximum probability')
    ax.set_zlabel('f(k,q,k)')
    ax.set_title(('f(k,q,k), 0<q<1 and 0<k<{:d} (l={:d})').
                 format(int(max(kVals)),length))

length = 128
nK = length*5
nQ = 150
q = np.linspace(0,1,nQ)

fig = pPlotUtil.figure()
ax1 = fig.add_subplot(1,2,1,projection='3d')
# plot for 'all k's', give a sense of 'global' decreasing
plotKMesh(length,np.linspace(0,length+1,nK),q,ax1)
ax2 = fig.add_subplot(1,2,2,projection='3d')
# plot for 'small k's', give a sense of detail decreasing
smallK = np.floor(np.log(length))
plotKMesh(length,np.linspace(0,smallK,nK),q,ax2)
pPlotUtil.savefig(fig,'q4')

