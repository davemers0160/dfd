import numpy as np
import math
#import cv2 as cv
from bokeh import events
# import bokeh
from bokeh.io import curdoc, output_file
from bokeh.models import ColumnDataSource, Spinner, Range1d, Slider, Legend, CustomJS, HoverTool, CheckboxGroup
from bokeh.plotting import figure, show, output_file
from bokeh.layouts import column, row, Spacer

from coc_calc import coc_calc


# f_num = 3.7
# f = 10
# d_o1 = 2
# d_o2 = 5
blur_img_w = 200
blur_img_h = 200
blur_alpha = np.full((blur_img_h, blur_img_w), 255, dtype=np.uint8)
blur_alpha[0:99, :] = 0
rs = np.arange(start=10, stop=1000, step=100)

limits = [0, 10000000]
step = 1000
spin_width = 110
# px_size = 0.00155

legend_label = ["fp 1", "fp 2", "CoC Difference"]

source = ColumnDataSource(data=dict(x=[], coc=[],  color=[], legend_label=[]))
stem_source = ColumnDataSource(data=dict(stem_x=[], stem_y0=[], stem_y1=[]))
#b_source = ColumnDataSource(data=dict(b1=[], b2=[]))
# source = ColumnDataSource(data=dict(x=[], coc=[], color=['blue', 'green', 'black'], legend_label=["fp 1", "fp 2", "CoC Difference"]))
# source = ColumnDataSource(data=dict(x=[], coc1=[], coc2=[], coc_diff=[]))

# setup the inputs
px_size = Spinner(title="pixel size (um)", low=0.001, high=10.0, step=0.001, value=3.45, width=spin_width)
f_num = Spinner(title="f number", low=0.01, high=100.0, step=0.01, value=2.35, width=spin_width)
f = Spinner(title="focal length (mm)", low=0.1, high=2500, step=0.1, value=200, width=spin_width)
#do_1 = Slider(title="focus point 1 (m):", start=1, end=20000, step=1, value=100, width=1400, callback_policy="mouseup", callback_throttle=50)
#do_2 = Slider(title="focus point 2 (m):", start=1, end=20000, step=1, value=10000, width=1400, callback_policy="mouseup")
do_1 = Slider(title="focus point 1 (m):", start=1, end=20000, step=1, value=100, width=1420)
do_2 = Slider(title="focus point 2 (m):", start=1, end=20000, step=1, value=10000, width=1420)
min_x_spin = Spinner(title="Min Range", low=limits[0], high=limits[1], step=1, value=0, width=spin_width)
max_x_spin = Spinner(title="Max Range", low=limits[0], high=limits[1], step=1, value=1000, width=spin_width)
min_y_spin = Spinner(title="Min Radius", low=limits[0], high=limits[1], step=1, value=0, width=spin_width)
max_y_spin = Spinner(title="Max Radius", low=limits[0], high=limits[1], step=1, value=40, width=spin_width)
div_spin = Spinner(title="Divisions", low=1, high=10000, step=1, value=10, width=spin_width)
hide_diff_cb = CheckboxGroup(labels=["Show CoC Diff"], active=[0], width=spin_width)
hide_dups_cb = CheckboxGroup(labels=["Show Duplicates"], active=[0], width=spin_width)

#fp1_b = figure(plot_height=200, plot_width=200, title="Focal Point 1", toolbar_location=None)
#fp1_b.image(image="b1", x=0, y=0, dw=200, dh=200, global_alpha=1.0, dilate=False, palette="Greys256", source=b_source)
#fp1_b.axis.visible = False
#fp1_b.grid.visible = False
#fp1_b.x_range.range_padding = 0
#fp1_b.y_range.range_padding = 0
#fp1_b.toolbar

ht = HoverTool(tooltips=[("Range: ", "$x"),  ("Blur Radius: ", "$y{0,0}")], point_policy="snap_to_data", mode="mouse")

coc_plot = figure(plot_height=550, plot_width=1300, title="Quantized Circles of Confusion")
l1 = coc_plot.multi_line(xs='x', ys='coc', source=source, line_width=2, color='color', legend_field='legend_label')
l2 = coc_plot.segment(x0='stem_x', y0='stem_y0', x1='stem_x', y1='stem_y1', source=stem_source, line_width=1, color='red')

coc_plot.xaxis.axis_label = "Range (m)"
coc_plot.yaxis.axis_label = "Pixel Radius"
coc_plot.axis.axis_label_text_font_style = "bold"
coc_plot.x_range = Range1d(start=min_x_spin.value, end=max_x_spin.value)
coc_plot.y_range = Range1d(start=min_y_spin.value, end=max_y_spin.value)
# coc_plot.legend[0].location = (900, 200))
# legend = Legend(items=[(("fp 1", "fp 2", "CoC Difference"), [l1])], location=(0, -60))
# coc_plot.add_layout(legend, 'right')
coc_plot.add_tools(ht)


# Custom JS code to update the plots
cb_dict = dict(source=source, stem_source=stem_source, coc_plot=coc_plot, px_size=px_size, f_num=f_num, f=f, do_1=do_1,
               do_2=do_2, min_x_spin=min_x_spin, max_x_spin=max_x_spin, min_y_spin=min_y_spin, max_y_spin=max_y_spin,
               div_spin=div_spin, limits=limits, step=step, hide_diff_cb=hide_diff_cb, hide_dups_cb=hide_dups_cb)
update_plot_callback = CustomJS(args=cb_dict, code="""
    var data = source.data;
    var stem_data = stem_source.data;
    
    data['x'] = [];
    data['coc'] = [];
    
    stem_data['stem_x'] = [];
    stem_data['stem_y0'] = [];
    stem_data['stem_y1'] = [];
    
    var d1 = do_1.value*1000;
    var d2 = do_2.value*1000;
    var px = px_size.value/1000;
    var start = Math.max(limits[0], step);
    var num = Math.floor((limits[1] - start) / step);
    
    var r = [];
    var coc1 = [];
    var coc2 = [];
    var coc_diff = [];
    
    var stem_x = [];
    var stem_y0 = [];
    var stem_y1 = [];

    var t1 = (d1*f.value*f.value)/(f_num.value*(d1-f.value));
    var t2 = (d2*f.value*f.value)/(f_num.value*(d2-f.value));
    
    var prev_coc = 0;
    var curr_coc = 0;
    
    console.log("test")
    
    for(var idx = 0; idx<num; idx++)
    {
        r.push((idx+1)*step)
        if (r[idx] < d1)
        {
            coc1.push(Math.ceil(t1*((1.0/r[idx]) - (1.0/d1))/px));
        }
        else if(r[idx] > d1)
        {
            coc1.push(Math.ceil(t1*((1.0/d1)-(1.0/r[idx]))/px));
        }
        else
        {
            coc1.push(0);
        }
        
        if (r[idx] < d2)
        {
            coc2.push(Math.ceil(t2*((1.0/r[idx]) - (1.0/d2))/px));
        }
        else if(r[idx] > d2)
        {
            coc2.push(Math.ceil(t2*((1.0/d2)-(1.0/r[idx]))/px));
        }
        else
        {
            coc2.push(0);
        }       
        
        coc_diff.push(Math.abs(coc1[idx]-coc2[idx]))
        r[idx] = r[idx]/1000;
        
        // get the ranges that correspond to the division settings
        if(r[idx]%div_spin.value == 0)
        {
            curr_coc = coc_diff[idx];
            if(curr_coc == prev_coc)
            {
                stem_x.push(r[idx]);
                stem_y0.push(0);
                stem_y1.push(max_y_spin.value);
            }
            
            prev_coc = curr_coc;
        }
       
    }
    
    data['x'].push(r);
    data['x'].push(r);
    data['x'].push(r);
    data['coc'].push(coc1);
    data['coc'].push(coc2);
    
    if(hide_diff_cb.active.length == 1)
    {
        data['coc'].push(coc_diff);
    }
    else
    {
        data['coc'].push([]);    
    }
    
    if(hide_dups_cb.active.length == 1)
    {
        stem_data['stem_x'] = stem_x;
        stem_data['stem_y0'] = stem_y0;
        stem_data['stem_y1'] = stem_y1;
    }
    else
    {
        stem_data['stem_x'] = [];
        stem_data['stem_y0'] = [];
        stem_data['stem_y1'] = [];    
    }
    
    coc_plot.x_range.start = min_x_spin.value;
    coc_plot.y_range.start = min_y_spin.value   
    coc_plot.x_range.end = max_x_spin.value;
    coc_plot.y_range.end = max_y_spin.value;
    
    source.change.emit();
    stem_source.change.emit();
    
    
""")


def update_plot(attr, old, new):

    r, coc1, coc_max1 = coc_calc(f_num.value, f.value, do_1.value * 1000, limits, step)
    r, coc2, coc_max2 = coc_calc(f_num.value, f.value, do_2.value * 1000, limits, step)

    r = r/1000
    px = px_size.value/1000

    q_coc1 = np.ceil(coc1/px)
    q_coc2 = np.ceil(coc2/px)

    q_coc_diff = abs(q_coc1 - q_coc2)

    r_slice = r[0:-1:div_spin.value]
    q_coc_slice = q_coc_diff[0:-1:div_spin.value]
    diff_match = np.hstack(([False], (q_coc_slice[0:-1] == q_coc_slice[1:])))

    source.data = dict(x=[r, r, r], coc=[q_coc1, q_coc2, q_coc_diff], color=['blue', 'green', 'black'], legend_label=["fp 1", "fp 2", "CoC Difference"])
    stem_source.data = dict(stem_x=r_slice[diff_match], stem_y0=np.zeros(r_slice[diff_match].shape[0]), stem_y1=np.full(r_slice[diff_match].shape[0], max_y_spin.value))
    # source.data = dict(x=[r], coc1=[q_coc1], coc2=[q_coc2], coc_diff=[q_coc_diff])
    #b_source.data = dict(b1=[blur_alpha], b2=[blur_alpha])

    coc_plot.x_range.start = min_x_spin.value
    coc_plot.y_range.start = min_y_spin.value

    coc_plot.x_range.end = max_x_spin.value
    coc_plot.y_range.end = max_y_spin.value


# setup the event callbacks for the plot
for w in [px_size, f_num, f, do_1, do_2, min_x_spin, max_x_spin, min_y_spin, max_y_spin, div_spin]:
    # w.on_change('value', update_plot)
    w.js_on_change('value', update_plot_callback)

for w in [hide_diff_cb, hide_dups_cb]:
    w.js_on_click(update_plot_callback)


update_plot(1, 1, 1)

# create the layout for the controls
range_input = column(min_x_spin, max_x_spin)
radius_input = column(min_y_spin, max_y_spin)

inputs = column(px_size, f_num, f, range_input, radius_input, div_spin, hide_diff_cb, hide_dups_cb)
layout = column(row(inputs, Spacer(width=20, height=20), coc_plot), do_1, do_2)

show(layout)

# doc = curdoc()
# doc.title = "Blur Calculator"
# doc.add_root(layout)

# output_file("d:/test.html", title='Bokeh Plot', mode='cdn', root_dir=None)
