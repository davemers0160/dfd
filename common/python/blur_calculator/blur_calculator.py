import numpy as np
import math
#import cv2 as cv
from bokeh import events
# import bokeh
from bokeh.io import curdoc, output_file
from bokeh.models import ColumnDataSource, Spinner, Range1d, Slider, Legend, CustomJS, HoverTool
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

limits = [0, 10000000]
step = 1000
spin_width = 100
# px_size = 0.00155

legend_label = ["fp 1", "fp 2", "CoC Difference"]

source = ColumnDataSource(data=dict(x=[], coc=[],  color=[], legend_label=[]))
b_source = ColumnDataSource(data=dict(b1=[], b2=[]))
# source = ColumnDataSource(data=dict(x=[], coc=[], color=['blue', 'green', 'black'], legend_label=["fp 1", "fp 2", "CoC Difference"]))
# source = ColumnDataSource(data=dict(x=[], coc1=[], coc2=[], coc_diff=[]))

# setup the inputs
px_size = Spinner(title="pixel size (um)", low=0.001, high=10.0, step=0.001, value=3.45, width=spin_width)
f_num = Spinner(title="f number", low=0.01, high=100.0, step=0.01, value=2.35, width=spin_width)
f = Spinner(title="focal length (mm)", low=0.1, high=2500, step=0.1, value=200, width=spin_width)
do_1 = Slider(title="focus point 1 (m):", start=1, end=20000, step=1, value=100, width=1400, callback_policy="mouseup", callback_throttle=50)
do_2 = Slider(title="focus point 2 (m):", start=1, end=20000, step=1, value=10000, width=1400, callback_policy="mouseup")
min_x_spin = Spinner(title="Min Range", low=limits[0], high=limits[1], step=1, value=0, width=spin_width)
max_x_spin = Spinner(title="Max Range", low=limits[0], high=limits[1], step=1, value=1000, width=spin_width)
min_y_spin = Spinner(title="Min Radius", low=limits[0], high=limits[1], step=1, value=0, width=spin_width)
max_y_spin = Spinner(title="Max Radius", low=limits[0], high=limits[1], step=1, value=40, width=spin_width)

fp1_b = figure(plot_height=200, plot_width=200, title="Focal Point 1", toolbar_location=None)
fp1_b.image(image="b1", x=0, y=0, dw=200, dh=200, global_alpha=1.0, dilate=False, palette="Greys256", source=b_source)
fp1_b.axis.visible = False
fp1_b.grid.visible = False
fp1_b.x_range.range_padding = 0
fp1_b.y_range.range_padding = 0
fp1_b.toolbar

ht = HoverTool(tooltips=[("Range: ", "$x"),  ("Blur Radius: ", "$y{0,0}")], point_policy="snap_to_data", mode="mouse")

coc_plot = figure(plot_height=550, plot_width=1300, title="Quantized Circles of Confusion")
l1=coc_plot.multi_line(xs='x', ys='coc', source=source, line_width=2, color='color', legend='legend_label')
# coc_plot.legend.title = "CoCs"
# coc_plot.line('x', 'coc1', source=source, line_width=2, color='blue', legend=legend_label[0])
# coc_plot.line('x', 'coc2', source=source, line_width=2, color='green', legend=legend_label[1])
# coc_plot.line('x', 'coc_diff', source=source, line_width=2, color='black', line_dash=(2,2), legend=legend_label[2])
coc_plot.xaxis.axis_label = "Range (m)"
coc_plot.yaxis.axis_label = "Pixel Radius"
coc_plot.axis.axis_label_text_font_style = "bold"
coc_plot.x_range = Range1d(start=min_x_spin.value, end=max_x_spin.value)
coc_plot.y_range = Range1d(start=min_y_spin.value, end=max_y_spin.value)
# coc_plot.legend[0].location = (900, 200))
# legend = Legend(items=[(("fp 1", "fp 2", "CoC Difference"), [l1])], location=(0, -60))
# coc_plot.add_layout(legend, 'right')
coc_plot.add_tools(ht)



b1_data = np.full((1, blur_img_h * blur_img_w), 255, dtype=np.uint8)
b1_data[:, 0:(blur_img_w * 99)] = 0
tmp_dict = dict(source=source, b_source=b_source, b1_d=b1_data)
tmp_cb = CustomJS(args=tmp_dict, code="""
    var data = source.data;
    
    //var b1_data = [];
    //Object.assign(b1_data, b1_d);
    var b1_data = JSON.parse(JSON.stringify(b1_d));
    
    //var b1_data = {...b1_d};
    //var b2_data = b_source.data['b2'];
    
    var stride = 200;
    var start = 10*stride;  
    var stop = 191*stride;
    
    x = Math.floor(cb_obj['x']);

    var b1 = data['coc'][0][x];
    var b2 = data['coc'][1][x];
    
    console.log(b1);
    console.log(b2);  
    
    for(var b=0; b<b1; b++)
    {
        for(var idx=start; idx<stop; idx++)
        {
            b1_data[0][idx] = 0.333333 * (b1_d[0][idx-stride] + b1_d[0][idx] + b1_d[0][idx+stride]);
            //for(var c=1; c<199; c++)
            //{
                //var p = b_source.data['b1']
                //b1_data[r][c] = Math.floor(0.111111*(b_source.data['b1'][r-1][c-1] + b_source.data['b1'][r-1][c] + b_source.data['b1'][r-1][c+1] + 
                //    b_source.data['b1'][r][c-1] + b_source.data['b1'][r][c] + b_source.data['b1'][r][c+1] + 
                //    b_source.data['b1'][r+1][c-1] + b_source.data['b1'][r+1][c] + b_source.data['b1'][r+1][c+1]));
            //}
        }
        
        for(var idx=start; idx<stop; idx++)
        {
            b1_data[0][idx] = Math.floor(0.333333 * (b1_d[0][idx-1] + b1_d[0][idx] + b1_d[0][idx+1]));
        }
    }
    
    
    b_source.data['b1'] = [];
    b_source.data['b1'].push(new Uint8Array(b1_data[0]));
    b_source.change.emit();
    
    """)


# Custom JS code to update the plots
cb_dict = dict(source=source, coc_plot=coc_plot, px_size=px_size, f_num=f_num, f=f, do_1=do_1, do_2=do_2,
               min_x_spin=min_x_spin, max_x_spin=max_x_spin, min_y_spin=min_y_spin, max_y_spin=max_y_spin,
               limits=limits, step=step)
update_plot_callback = CustomJS(args=cb_dict, code="""
    var data = source.data;
    data['x'] = [];
    data['coc'] = [];
    
    var d1 = do_1.value*1000;
    var d2 = do_2.value*1000;
    var px = px_size.value/1000;
    var start = Math.max(limits[0], step);
    var num = Math.floor((limits[1] - start) / step);
    
    var r = [];
    var coc1 = [];
    var coc2 = [];
    var coc_diff = [];
    
    var t1 = (d1*f.value*f.value)/(f_num.value*(d1-f.value));
    var t2 = (d2*f.value*f.value)/(f_num.value*(d2-f.value));
    
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
       
    }
    
    data['x'].push(r);
    data['x'].push(r);
    data['x'].push(r);
    data['coc'].push(coc1);
    data['coc'].push(coc2);
    data['coc'].push(coc_diff);
    
    coc_plot.x_range.start = min_x_spin.value;
    coc_plot.y_range.start = min_y_spin.value   
    coc_plot.x_range.end = max_x_spin.value;
    coc_plot.y_range.end = max_y_spin.value
    source.change.emit();
""")


def update_plot(attr, old, new):

    r, coc1, coc_max1 = coc_calc(f_num.value, f.value, do_1.value * 1000, limits, step)
    r, coc2, coc_max2 = coc_calc(f_num.value, f.value, do_2.value * 1000, limits, step)

    r = r/1000
    px = px_size.value/1000

    q_coc1 = np.ceil(coc1/px)
    q_coc2 = np.ceil(coc2/px)

    q_coc_diff = abs(q_coc1 - q_coc2)

    source.data = dict(x=[r, r, r], coc=[q_coc1, q_coc2, q_coc_diff], color=['blue', 'green', 'black'], legend_label=["fp 1", "fp 2", "CoC Difference"])
    # source.data = dict(x=[r], coc1=[q_coc1], coc2=[q_coc2], coc_diff=[q_coc_diff])
    b_source.data = dict(b1=[blur_alpha], b2=[blur_alpha])

    coc_plot.x_range.start = min_x_spin.value
    coc_plot.y_range.start = min_y_spin.value

    coc_plot.x_range.end = max_x_spin.value
    coc_plot.y_range.end = max_y_spin.value


for w in [px_size, f_num, f, do_1, do_2, min_x_spin, max_x_spin, min_y_spin, max_y_spin]:
    # w.on_change('value', update_plot)
    w.js_on_change('value', update_plot_callback)

# coc_plot.js_on_event(events.Tap, tmp_cb)
# coc_plot.js_on_event(events.Tap, blur_cb)

update_plot(1, 1, 1)

# layout = column([row([max_x_spin,max_y_spin]), blur_plot])

range_input = column(min_x_spin, max_x_spin)
radius_input = column(min_y_spin, max_y_spin)
#range_pn = Panel(child=range_input, title="Range")

inputs = column(px_size, f_num, f, range_input, radius_input)
#blur_imgs = column(row([Spacer(width=20, height=20), fp1_b]), row([Spacer(width=20, height=20), fp1_b]))
layout = column(row(inputs, Spacer(width=20, height=20), coc_plot), do_1, do_2)

show(layout)

# doc = curdoc()
# doc.title = "Blur Calculator"
# doc.add_root(layout)

# output_file("d:/test.html", title='Bokeh Plot', mode='cdn', root_dir=None)
