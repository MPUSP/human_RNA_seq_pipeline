# Interactive Scatter plot from DESeq2 results
# 21.01.20
# Tim Sullivan

import sys

import pandas as pd
import numpy as np 

from bokeh.io import output_file, show, save
from bokeh.models import ColumnDataSource, HoverTool, LinearColorMapper, CDSView, BooleanFilter, Panel, Tabs, CustomJS, Slider, RadioGroup
from bokeh.models.widgets import Div
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.plotting import figure
from bokeh.layouts import column, row

#Data loading section
import os
print (os.getcwd())



deseq_df = pd.read_csv(sys.argv[1], sep=',')
comparison_name = sys.argv[2]


# spydf.sort_values(['padj'], inplace=True, axis=0)

deseq_df = deseq_df[deseq_df.padj < .9]

deseq_df = deseq_df.dropna(axis=0)
deseq_df.reset_index(inplace=True)

def calc_colors(pvals, ef, limit, ef_limit):
    color_list = ['gray']*len(pvals)
    
    c = pvals<limit

    for i in range(len(c)):
        if abs(ef[i]) < ef_limit:
            continue
        if c[i]:
            color_list[i] = 'red'
    return color_list
    

source = ColumnDataSource(data=dict(l2fc = deseq_df['log2FoldChange'],
                                    desc = deseq_df['GENE'],
                                    pval = deseq_df['padj']+10**-300,
                                    minuspval = -1*np.log10(deseq_df['padj']+10**-300),
                                    mean = deseq_df['baseMean'],
                                    l2reads = np.log2(deseq_df['baseMean']),
                                    pointsize = 6+(2*np.log10(deseq_df['baseMean']+.1)),
                                    pointcolor = calc_colors(deseq_df['padj'], deseq_df['log2FoldChange'], .05, 1.0),
                                    ))

    

#Scatterplot Section
hover_tooltips = [
    ('Gene', '@desc'),
    ('P-value', '@pval'),
    ('Log2FC', '@l2fc'),
    ('Mean Read Counts', '@mean'),
]


p = figure(plot_width=800, plot_height=600,
           tools="hover,  box_zoom, pan, reset, save",
           active_drag = "box_zoom",
           tooltips = hover_tooltips,
           x_range=(-15,15), y_range=(0,max(-1*np.log10(deseq_df['padj']+10**-300))+5),
           x_axis_label = "Log2 Fold Change",
           y_axis_label = "-log10(pvalue)",
            )

p.circle('l2fc', 'minuspval', size='pointsize', fill_alpha=.6, line_alpha=.4, source=source, fill_color='pointcolor', line_color="gray")
# p.line([-15,15],[-15,15], color='pink', alpha=.6)
# p.line([-15,15],[-14,16], color='red', alpha=0)
# p.line([-15,15],[-16,14], color='red', alpha=0)
p.xaxis.axis_label_text_font_size = "18pt"
p.yaxis.axis_label_text_font_size = "18pt"
p.xaxis.axis_label_text_font_style = 'normal'
p.yaxis.axis_label_text_font_style = 'normal'
p.title.text_font_size='24pt'
p.toolbar.logo=None





pval_slider = Slider(title="-log10(p-value) coloring cutoff", start=0, end=10, value=2, step=1)
# color_buttons = RadioGroup(labels=["One Color", "Three Colors"], active=1)
ef_slider = Slider(title="log2FC effect size difference", start=0, end=10, value=1, step=1)


pval_callback = CustomJS(args=dict(source=source, pval=pval_slider, ef=ef_slider), code="""                   
    const data = source.data;
    const pval_limit = 10**(-1*pval.value);
    const ef_limit = ef.value;
    const pointcolor = data['pointcolor'];
    
    const p = data['pval'];    
    const ef_data = data['l2fc'];

    console.log(ef_limit)
    console.log(ef_data)
    console.log(pointcolor)
    
    for (var i=0, len=pointcolor.length; i < len; i++ ){
            pointcolor[i] = 'gray'   
            if (Math.abs(ef_data[i])< ef_limit){  continue; }
            
            if ( p[i]<pval_limit){pointcolor[i] = 'red';continue;}
 
            }
  
    source.change.emit();
""")



pval_slider.js_on_change('value', pval_callback)
ef_slider.js_on_change('value', pval_callback)
# color_buttons.js_on_change('active', pval_callback)


controls_div = Div(text="""
        <h2>Controls</h2>
""",
width=900, height=45)

# plot_div = Div(text="""
#         <h2>Log2FC Plot</h2>
# """,
# width=900, height=55)
# layout = column(
#     controls_div,
#     ef_slider,
#     pval_slider,
#     color_buttons,
#     plot_div,
#     p
# )
plot_div = Div(text="""
        <h2>Volcano Plot</h2>
""",
width=900, height=55)
layout = (column(    controls_div,
    ef_slider,
    pval_slider,plot_div,p))

tab1 = Panel(child=layout, title='Volcano Plot')




p = figure(plot_width=800, plot_height=600,
           tools="hover,  box_zoom, pan, reset, save",
           active_drag = "box_zoom",
           tooltips = hover_tooltips,
           x_range=(3,18), y_range=(-8,15),
           x_axis_label = "Log2 Read Count",
           y_axis_label = "log2 Fold Change",
            )

p.circle('l2reads', 'l2fc', size=8, fill_alpha=.6, line_alpha=.4, source=source, fill_color='pointcolor', line_color="gray")
# p.line([-15,15],[-15,15], color='pink', alpha=.6)
# p.line([-15,15],[-14,16], color='red', alpha=0)
# p.line([-15,15],[-16,14], color='red', alpha=0)
p.xaxis.axis_label_text_font_size = "18pt"
p.yaxis.axis_label_text_font_size = "18pt"
p.xaxis.axis_label_text_font_style = 'normal'
p.yaxis.axis_label_text_font_style = 'normal'
p.title.text_font_size='24pt'
p.toolbar.logo=None




plot_div = Div(text="""
        <h2>MA Plot</h2>
""",
width=900, height=55)
layout = (column(    controls_div,
    ef_slider,
    pval_slider,plot_div,p))

tab2 = Panel(child=layout, title='MA Plot')



#Readme Section
p = Div(text="""
<h1> Interactive Volcano and MA plots</h1>

<h4> The source of these data are the samples from SARS-CoV-2 vs mock infected cell lines, from external projects PRJNA615032 or PRJNA625518</h4>

Here are a <a href="https://en.wikipedia.org/wiki/Volcano_plot_(statistics)">Volcano plot</a> and an <a href="https://en.wikipedia.org/wiki/MA_plot">MA plot</a> used to explore relationships between p-values, fold changes, and expression levels of genes in differentially expression experiments.
<br><br>
Hovering over each point will show the multiple-hypothesis corrected p-values for each test, the log2 fold change, and the mean expression level across the replicates in each treatment. 
<br><br>
(<strong>p-values more significant than 10^-300</strong> have been set to 10^-300)<br>
<br><br>
The buttons on the right toolbar of the plot can be used to zoom by selecting a group of points, or pan the image.
<br><br>
<font color="red"><strong>Warning!</strong></font> Hovering over too many points at once will cause the page to freeze!  If this happens, refreshing to page will fix it.

</p>
<p>
<br><small> Generated by Tim Sullivan on 29.09.2020 </small>
</p>

""",
width=900, height=700)
tab3 = Panel(child=p, title='Readme')


#Raw Data Section
data_columns = [
        TableColumn(field="desc", title="Gene ID"),
        TableColumn(field="l2fc", title="log2FC",),
        TableColumn(field="pval", title="adj. P-value"),
        TableColumn(field="mean", title="Mean Reads"),
    ]
data_table = DataTable(source=source, columns=data_columns, width=900, height=700)
tab4 = Panel(child=data_table, title='Raw Data')


#Tab convolution
tabs = Tabs(tabs=[tab1, tab2, tab3, tab4])
#Output
output_file('interactive_volcano_ma_'+comparison_name+'.html')
print('interactive_volcano_ma_'+comparison_name+'.html')
save(tabs)




