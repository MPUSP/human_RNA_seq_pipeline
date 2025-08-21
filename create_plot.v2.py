#!./venv/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster


class defaults:
    
    class input_data:
        mv_data_file = ['result_tables/H4_UvsM_DEGs.csv', 'result_tables/MVs_vs_Unt.csv']
        spy_data_file = ['result_tables/H4_UvsS_DEGs.csv', 'result_tables/Spy_vs_Unt.csv']
    
    pval_threshold = 0.01
    fc_threshold = 2
    fc_ratio_threshold = 4#in lin space
    
    selected_genes = ['GBP5', 'IDO1', 'CXCL9', 'CXCL11', 'IFNB1', 'CXCL10', 'IL36G', 'IL12B', 'CCL20', 'CSF3', 'CCL4', 'IL6', 'IL1B', 'CCL3', 'TNF']#sullivan, BLA cells
    # selected_genes = ['CCR2', 'PTGFRN', 'FOSB', "CXCL10", 'IFNG', 'IFNB1', 'TNF', 'IL1A', 'IL12B', 'CSF3', 'IL6', 'IL36G', 'CCL20', 'CCL4', 'IL1B', 'IL19', 'WNT5B', 'CASP5', 'WNT5A', 'STEAP4']# eric
    
    class color:
        base = (120/255, 132/255, 131/255,.2)
        p_mv = (252/255,147/255,0/255,.3)
        p_spy = (36/255,0/255,156/255,.3)
        p_shared = (0/255,187/255,250/255,.3)
        grid = (.9,.9,.9)
    
    axis = [-11, 16]
    dist_threshold = 1.25#for annotation of genes
    anno_offset = 15
    anno_size = 6
    plot_area_lines = False
    
        

def read_data(datei): #{geneID: (baseMean:float, log2FoldChange:float, lfcSE:float, pvalue:float, padj:float, Gene:str)}
    dataset = {}
    with open(datei, 'r') as infile:
        head = False
        for line in infile:
            if head:
                line = line.strip().split(',')
                try:
                    dataset[line[1]] = [float(line[v]) for v in range(2, 8)]+[line[8]]
                except:
                    pass
            else:
                head = True
    return dataset

def read_meta(datei): #{sampleID: (treatment:str, run:str)}
    dataset = {}
    with open(datei, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            dataset[line[0]] = (line[1], line[2])
    return dataset

def create_scatter_plot(mv, spy, prefix=''):
    #create scatter plot
    fig, ax = plt.subplots(figsize=(5,5))
    gene_list = set(spy.keys()) & set(mv.keys())
    
    x_vals_base = []
    y_vals_base = []
    x_vals_p_mv = []
    y_vals_p_mv = []
    x_vals_p_spy = []
    y_vals_p_spy = []
    x_vals_p_shared = []
    y_vals_p_shared = []
    
    x_vals_selected = []
    y_vals_selected = []
    names_selected = []
    
    x_size_base = []
    y_size_base = []
    x_size_mv = []
    y_size_mv = []
    x_size_spy = []
    y_size_spy = []
    x_size_shared = []
    y_size_shared = []
    
    size_column = 0
    
    for gene in gene_list:
        if 1/defaults.fc_ratio_threshold < 2**mv[gene][1]/2**spy[gene][1] < defaults.fc_ratio_threshold:
            x_vals_base.append(spy[gene][1])
            y_vals_base.append(mv[gene][1])
            x_size_base.append(spy[gene][size_column])
            y_size_base.append(mv[gene][size_column])
        else:
            if mv[gene][4] < defaults.pval_threshold and spy[gene][4] < defaults.pval_threshold and abs(mv[gene][1]) > defaults.fc_threshold and abs(spy[gene][1]) > defaults.fc_threshold:
                x_vals_p_shared.append(spy[gene][1])
                y_vals_p_shared.append(mv[gene][1])
                x_size_shared.append(spy[gene][size_column])
                y_size_shared.append(mv[gene][size_column])
            elif mv[gene][4] < defaults.pval_threshold and abs(mv[gene][1]) > defaults.fc_threshold:
                x_vals_p_mv.append(spy[gene][1])
                y_vals_p_mv.append(mv[gene][1])
                x_size_mv.append(spy[gene][size_column])
                y_size_mv.append(mv[gene][size_column])
            elif spy[gene][4] < defaults.pval_threshold and abs(spy[gene][1]) > defaults.fc_threshold:
                x_vals_p_spy.append(spy[gene][1])
                y_vals_p_spy.append(mv[gene][1])
                x_size_spy.append(spy[gene][size_column])
                y_size_spy.append(mv[gene][size_column])
            else:
                x_vals_base.append(spy[gene][1])
                y_vals_base.append(mv[gene][1])
                x_size_base.append(spy[gene][size_column])
                y_size_base.append(mv[gene][size_column])
                
        if spy[gene][6] in defaults.selected_genes:
            x_vals_selected.append(spy[gene][1])
            y_vals_selected.append(mv[gene][1])
            names_selected.append(spy[gene][6])
            
    
    ax.grid(True, color=defaults.color.grid, zorder=0)
    
    ax.scatter(x_vals_base, y_vals_base, color=defaults.color.base, edgecolor=None, lw=0, s=5, zorder=100)
    ax.scatter(x_vals_p_mv, y_vals_p_mv, color=defaults.color.p_mv, edgecolor=None, lw=0, s=5, zorder=101)
    ax.scatter(x_vals_p_spy, y_vals_p_spy, color=defaults.color.p_spy, edgecolor=None, lw=0, s=5, zorder=102)
    ax.scatter(x_vals_p_shared, y_vals_p_shared, color=defaults.color.p_shared, edgecolor=None, lw=0, s=5, zorder=103)
    
    #adding lines
    ##x==y
    ax.plot([defaults.axis[0], defaults.axis[1]], [defaults.axis[0], defaults.axis[1]], color=defaults.color.grid, lw=.5, zorder=200)
    
    if defaults.plot_area_lines:
        ##lowerleft
        ax.plot([defaults.axis[0], -defaults.fc_threshold], [-defaults.fc_threshold, -defaults.fc_threshold], color='black', lw=1, ls='--', zorder=201)
        ax.plot([-defaults.fc_threshold, -defaults.fc_threshold], [defaults.axis[0], -defaults.fc_threshold], color='black', lw=1, ls='--', zorder=202)
        
        ##upperright
        ax.plot([5, defaults.axis[1]], [5,5], color='black', lw=1, ls='--', zorder=203)
        ax.plot([5, 5], [5, defaults.axis[1]], color='black', lw=1, ls='--', zorder=204)
        
        ##middle
        ax.plot([defaults.axis[0], defaults.fc_threshold-math.log(defaults.fc_ratio_threshold, 2)], [defaults.fc_threshold, defaults.fc_threshold], color='black', lw=1, ls='--', zorder=205)
        ax.plot([defaults.fc_threshold, defaults.fc_threshold], [defaults.axis[0], defaults.fc_threshold-math.log(defaults.fc_ratio_threshold, 2)], color='black', lw=1, ls='--', zorder=206)
        
        ##diagonals
        ax.plot([defaults.fc_threshold-math.log(defaults.fc_ratio_threshold, 2), 5], [defaults.fc_threshold, 5+ math.log(defaults.fc_ratio_threshold, 2)], color='black', lw=1, ls='--', zorder=207)
        ax.plot([defaults.fc_threshold, 5+ math.log(defaults.fc_ratio_threshold, 2)], [defaults.fc_threshold-math.log(defaults.fc_ratio_threshold, 2), 5], color='black', lw=1, ls='--', zorder=208)
    
        #add quadrand text
        for e, coord in enumerate([(defaults.axis[0]+1, defaults.axis[1]-1), (defaults.axis[1]-1, defaults.axis[1]-1), (defaults.axis[1]-1, defaults.axis[0]+1), (defaults.axis[0]+1, defaults.axis[0]+1)]):
            ax.text(coord[0], coord[1], "$\mathbf{"+f"{str(e+1)}"+"}$", ha='center', va='center', fontsize=10, color='black', zorder=300)
    
    
    ##annotate selected genes
    plot_groups = [[], []]
    for x, y, name in zip(x_vals_selected, y_vals_selected, names_selected):
        if y>x: plot_groups[0].append((x,y,name))
        else: plot_groups[1].append((x,y,name))
    
    lookup = {n:(x,y) for x,y,n in zip(x_vals_selected, y_vals_selected, names_selected)}
    
    for sub in plot_groups:
        data_dist = pdist(np.array([(x,y) for x,y,n in sub]))
        data_link = linkage(data_dist)
        clusters = fcluster(linkage(data_dist), defaults.dist_threshold, criterion='distance')
        groups = {}
        for cluster, name in zip(clusters, [n for x,y,n in sub]):
            if cluster in groups: groups[cluster].append(name)
            else: groups[cluster] = [name]

        for c in groups.keys():
            if len(groups[c]) == 1:
                group_center = (lookup[groups[c][0]][0] + lookup[groups[c][0]][1])/2
                
            elif len(groups[c]) > 1:
                group_center = (np.mean([lookup[n][0] for n in groups[c]]), np.mean([lookup[n][1] for n in groups[c]]))
                group_center = (group_center[0] + group_center[1])/2
                
            for name in groups[c]:
                ha_text = 'right'
                va_text = 'baseline'
                offset = (-defaults.anno_offset, defaults.anno_offset)
                if lookup[name][0] > lookup[name][1]:
                    ha_text = 'left'
                    va_text = 'top'
                    offset = (defaults.anno_offset, -defaults.anno_offset)
                offset = radial_offset(lookup[name][0], lookup[name][1], group_center, group_center, defaults.anno_offset)
                plt.annotate(name, (lookup[name][0], lookup[name][1]), textcoords='offset points', xytext=offset, ha=ha_text, va=va_text, fontsize=defaults.anno_size, color='black', arrowprops=dict(arrowstyle="-", color='black', lw=.5), zorder=301, bbox=dict(boxstyle="round, pad=0.05", edgecolor=(1,1,1,0), facecolor=(1,1,1,.5)))
    
    # for x, y, name in zip(x_vals_selected, y_vals_selected, names_selected):
    #     ha_text = 'right'
    #     va_text = 'baseline'
    #     offset = (-defaults.anno_offset, defaults.anno_offset)
    #     if x > y:
    #         ha_text = 'left'
    #         va_text = 'top'
    #         offset = (defaults.anno_offset, -defaults.anno_offset)
    #     plt.annotate(name, (x, y), textcoords='offset points', xytext=offset, ha=ha_text, va=va_text, fontsize=defaults.anno_size, color='black', arrowprops=dict(arrowstyle="-", color='black', lw=.5), zorder=301, bbox=dict(boxstyle="round, pad=0.05", edgecolor=(1,1,1,0), facecolor=(1,1,1,.5)))
        
    
    #formatting figure
    ax.set_xlabel('$\mathit{Spy}$ vs Unt Log$_{2}$FC')
    ax.set_ylabel('EVs vs Unt Log$_{2}$FC')
    
    ax.set_xlim(defaults.axis[0], defaults.axis[1])
    ax.set_ylim(defaults.axis[0], defaults.axis[1])
    
    plt.tight_layout()
    
    for typ in ['svg', 'png']:
        plt.savefig(f'{prefix}scatter_plot.{typ}', dpi=300, format=typ)

def calc_dist_to_center(x, y, cx, cy):
    return np.sqrt((x-cx)**2+(y-cy)**2)

def radial_offset(x_point, y_point, center_x, center_y, distance):
        angle = np.arctan2(y_point - center_y, x_point - center_x)
        offset_x = np.cos(angle) * distance
        offset_y = np.sin(angle) * distance
        return offset_x, offset_y

def main():
    for i in range(2):
        prefix = defaults.input_data.mv_data_file[i].split('/')[-1].split('.')[0]+'_'
        
        #read data
        mv = read_data(defaults.input_data.mv_data_file[i])
        spy = read_data(defaults.input_data.spy_data_file[i])
            
        #create scatter plot
        create_scatter_plot(mv, spy)

main()
print('Done')