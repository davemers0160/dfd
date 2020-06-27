#
import platform
import os
import math
#import ctypes as ct
#from cffi import FFI
import numpy as np
import cv2 as cv
import bokeh
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, show
from bokeh.layouts import column, row, Spacer

script_path = os.path.realpath(__file__)

# set up some global variables that will be used throughout the code
image_name = os.path.dirname(script_path) + "/image_f1_000_100.png"
depth_name = os.path.dirname(script_path) + "/dm_000.png"


color_img = cv.cvtColor(cv.imread(image_name), cv.COLOR_BGR2RGBA)
depth_img = cv.cvtColor(cv.imread(depth_name), cv.COLOR_BGR2GRAY)

source = ColumnDataSource(data=dict(color_img=[color_img], depth_map=[depth_img]))

TOOLTIPS = [("value", "@depth_map")]

p1 = figure(plot_height=500, plot_width=500, title="Input", toolbar_location="right", tooltips=TOOLTIPS)
p1.image_rgba(image="color_img", x=0, y=0, dw=400, dh=400, source=source)
p1.axis.visible = False
p1.grid.visible = False
p1.x_range.range_padding = 0
p1.y_range.range_padding = 0


layout = column(p1)

show(layout)

