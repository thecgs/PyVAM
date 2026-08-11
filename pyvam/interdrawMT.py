#!/usr/bin/env python
# coding: utf-8

import math
import re
import logging
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd

from .parserGB import get_features
from .config import MTColors_legends, MTColors, FullName2AbbrName

class InteractiveMTVisualizer:
    def __init__(self):
        #self.colors = MTColors_legends
        self.default_layout = {'hovermode': 'closest',
                               'plot_bgcolor': 'white',
                               'paper_bgcolor': 'white',
                              }
        
        self.nav_config = {"editable": False, 
                           'displaylogo': False,
                           'displayModeBar': True,
                           'modeBarButtons': [['zoomIn2d', 'zoomOut2d', 'zoom2d', 'autoScale2d', 
                                               'resetScale2d', 'pan2d', 'select2d','lasso2d'],
                                              ['drawcircle', 'drawrect', 'drawline', 'drawopenpath',
                                               'drawclosedpath', 'eraseshape', 'toImage']
                                             ],
                           'toImageButtonOptions': {'format': 'svg'}
                          }
        
    def _convert_color(self, color_str):
        if not color_str or not isinstance(color_str, str):
            return color_str
            
        if not color_str.startswith('#'):
            return color_str
            
        if len(color_str) == 9:
            r = int(color_str[1:3], 16)
            g = int(color_str[3:5], 16)
            b = int(color_str[5:7], 16)
            a = int(color_str[7:9], 16) / 255.0
            return f'rgba({r}, {g}, {b}, {a:.2f})'
        elif len(color_str) == 7:
            return color_str
        elif len(color_str) > 9:
            return color_str[:7]
        else:
            return color_str
    
    def _get_box_param(self, gene_name):
        head_length=0.0025
        if gene_name in ["F","V","L","I","Q","M","W","A","N","C","Y","S","D","K","G","R","H","S","L","E","T","P"]:
            tail_length = 0.005
            text_offset = tail_length/2.5
            
        elif gene_name in ["ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"]:
            tail_length = 0.008
            if gene_name=="ND4L":
                text_offset = tail_length/8
            else:
                text_offset = tail_length/6
        elif gene_name in ["12S", "16S"]:
            tail_length = 0.005
            text_offset = tail_length/6
        elif gene_name in ['COX1','COX2','COX3']:
            tail_length = 0.008
            text_offset = tail_length/8
        elif gene_name in ["ATP6", "ATP8"]:
            tail_length = 0.008
            text_offset = tail_length/6
        elif gene_name == "NCR":
            tail_length = 0.007
            text_offset = tail_length/6
        elif gene_name == "Cytb":
            tail_length = 0.008
            text_offset = tail_length/5
        elif gene_name == 'Gap':
            tail_length = 0.005
            text_offset = tail_length/6
        else:
            tail_length = 0.015
            text_offset = tail_length/4
        box_width = head_length + tail_length
        return head_length, tail_length, box_width, text_offset

    def _get_colors_group(self, colors):
    
    
        res = {}
        if colors.upper() == "MITOZ":
            for gene in MTColors[colors.upper()]:
                if gene in ['ND1', 'ND2', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'COX1', 'COX2', 'COX3', 'ATPase6', 'ATPase8', 'Cytb']:
                    res[gene] = {'group':'Protein codon genes'}
                elif gene in ['tRNA-His', 'tRNA-Pro', 'tRNA-Thr', 'tRNA-Trp', 'tRNA-Met',
                            'tRNA-Asp', 'tRNA-Ala', 'tRNA-Gln', 'tRNA-Ile', 'tRNA-Arg',
                            'tRNA-Tyr', 'tRNA-Phe', 'tRNA-Lys', 'tRNA-Gly', 'tRNA-Asn',
                            'tRNA-Leu', 'tRNA-Glu', 'tRNA-Val', 'tRNA-Cys', 'tRNA-Ser']:
                    res[gene] = {'group':'transfer RNA genes'}
                elif gene in ['12S rRNA', '16S rRNA']:
                    res[gene] = {'group':'ribosomal RNA genes'}
                elif gene in ['D-loop']:
                    res[gene] = {'group':'D-loop (Non-coding region)'}
                else:
                    res[gene] = {'group':'Other genes'}
                    
        elif colors.upper() in ["CHEN", "OGDRAW", "MITOFISH", "MITOFISH1", "TAN", "CHLOROPLOT", "GREY", "GGGENES"]:
            for gene in MTColors[colors.upper()]:
                if gene in ['ND1', 'ND2', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6']:
                    res[gene] = {'group':'Complex I (NADH dehydrogenase)'}
                elif gene in ['COX1', 'COX2', 'COX3']:
                    res[gene] = {'group':'Complex IV (Cytochrome c oxidase)'}
                elif gene in ['ATPase6', 'ATPase8']:
                    res[gene] = {'group':'ATP synthase'}
                elif gene in ['Cytb']:
                    res[gene] = {'group':'Cytochrome b'}
                elif gene in ['tRNA-His', 'tRNA-Pro', 'tRNA-Thr', 'tRNA-Trp', 'tRNA-Met',
                            'tRNA-Asp', 'tRNA-Ala', 'tRNA-Gln', 'tRNA-Ile', 'tRNA-Arg',
                            'tRNA-Tyr', 'tRNA-Phe', 'tRNA-Lys', 'tRNA-Gly', 'tRNA-Asn',
                            'tRNA-Leu', 'tRNA-Glu', 'tRNA-Val', 'tRNA-Cys', 'tRNA-Ser']:
                    res[gene] = {'group':'transfer RNA'}
                elif gene in ['12S rRNA', '16S rRNA']:
                    res[gene] = {'group':'ribosomal RNA'}
                elif gene in ['D-loop']:
                    res[gene] = {'group':'D-loop (Non-coding region)'}
                else:
                    res[gene] = {'group':'Other genes'}
        else:
            if colors.upper() == "IGV":
                for gene in MTColors[colors.upper()]:
                    res[gene] = {'group':gene}
            else:
                for gene in colors:
                    res[gene] = {'group':gene}

        res['Gap'] =  {'group':"Gap"}
        return res
        
    def remove_join(self, features):
        tmp = []
        res = []
        for feature in features:
            if feature.join!=None:
                if feature.join not in tmp:
                    res.append(feature)
                    tmp.append(feature.join)
            else:
                res.append(feature)
        return res
    
    def draw_linear_MT_nonproportional_plotly(self,
                                              files,
                                              output=None,
                                              abbr=False,
                                              isfilename2species=False,
                                              colors="mitofish",
                                              gene_label_size=9,
                                              gene_label_color='black',
                                              show_legend=True,
                                              start=None,
                                              species_label_size=12, 
                                              species_label_color='black',
                                              add_id=False,
                                              editable=False,
                                              force_reoriented=False,
                                              remove_NCR = False,
                                              ):
        """
        Descripton:
            The order of mitochondrial genes is not depicted in proportion to their gene size.
            Interactive chart developed using Plotly.
        
        Parameters：
            files: {str, list, tuple} one or more genbankfile or NCBI accession ID.
            output: {str} a path of fig save, save as HTML file.
            abbr: {bool} whether to abbreviate species names.
            isfilename2species: {bool} whether filename convert to species.
            colors: {str, dict} themes such as, Chen, Tan, ogdraw, mitofish,
                                mitofish1, mitoz,  gggenes, chloroplot, grey, igv.
            gene_label_size: {int} gene label size.
            gene_label_color: {str} gene label color.
            show_legend: {bool} show fig legend.
            start: {None, str} initial feature, such as, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
                         COX1, COX2, COX3, ATPase6, ATPase8, Cytb, tRNA-His, tRNA-Pro,
                         tRNA-Thr, tRNA-Trp, tRNA-Met, tRNA-Asp, tRNA-Ala, tRNA-Gln,
                         tRNA-Ile, tRNA-Arg, tRNA-Tyr, tRNA-Phe, tRNA-Lys, tRNA-Gly,
                         tRNA-Asn, tRNA-Leu, tRNA-Glu, tRNA-Val, tRNA-Cys, tRNA-Ser,
                         12S rRNA, 16S rRNA, D-loop. default=None.
            species_label_size: {int} species name size.
            species_label_color: {str} species name color.
            add_id: {bool} Species add to accession id from NCBI.
            editable: {bool} Make HTML files editable.
            force_reoriented: {bool} force-reoriendted linear mtgenome.
            remove_NCR {bool}: Do not display D-loop; this may be useful for comparing a mix of genbank with and without NCR annotation.
        """
        color2groups = self._get_colors_group(colors)
        tmp_groups = set()
        
        if editable:
            self.nav_config["editable"] = True
        
        if not isinstance(files, (list, tuple)):
            files = [files]
        
        fig = go.Figure()
        
        genomes = []
        for file in reversed(files):
            features = get_features(file, abbr=abbr, isfilename2species=isfilename2species, colors=colors, start=start, force_reoriented=force_reoriented)
            if remove_NCR:
                features = [feature for feature in features if feature.name != "D-loop"]
            features = self.remove_join(features)
            genomes.append(features)
        
        y_positions = []
        species_labels = []


        _staus_brake = False
        for y, features in enumerate(genomes):
            y_pos = y + 0.5
            y_positions.append(y_pos)
            
            species_label = features[0].name + (f" ({features[0].accession})" if add_id else "")
            species_labels.append(species_label)
            
            if features[1].location.strand == -1:
                x_position = self._get_box_param(features[1].name)[2]+0
            else:
                x_position=0
        
            for i, feature in enumerate(features[1:]):
                start_pos = int(feature.location.start) + 1
                end_pos = int(feature.location.end)
                length = end_pos - start_pos + 1
                
                if feature.join == None:
                    hover_text = f"""Gene Name: {feature.name}<br>Position: {start_pos:}-{end_pos}{'(+) ' if feature.location.strand == 1 else '(-) '}"""
                else:
                    hover_text = f"""Gene Name: {feature.name}<br>Position: {str(feature.join)}"""
                    
                head_length, tail_length, box_width, text_offset = self._get_box_param(FullName2AbbrName.get(feature.name, feature.name))
                fill_color = self._convert_color(feature.color)
                height = 0.5
                
                #if color2groups[feature.name].get('group', "Other") not in tmp_groups:
                #    showlegend = True
                #else:
                #    showlegend = False
                
                group_name = color2groups[feature.name]['group'] if feature.name in color2groups else color2groups["Other genes"]['group']
                if group_name not in tmp_groups:
                    showlegend = True
                else:
                    showlegend = False
                
                if feature.location.strand == 1:
                    fig.add_trace(go.Scatter(
                        x=[x_position, 
                           x_position + tail_length,
                           x_position + tail_length + head_length, 
                           x_position + tail_length,
                           x_position,
                           x_position,],
                        y=[y_pos - height/2,
                           y_pos - height/2,
                           y_pos,
                           y_pos + height/2,
                           y_pos + height/2,
                           y_pos - height/2],
                        mode='lines',
                        fill="toself",
                        fillcolor=fill_color,
                        line=dict(color='black', width=1),
                        name=group_name,
                        legendgroup=group_name,
                        hoverinfo='skip',
                        #hoverinfo='text',
                        #hoveron='points+fills',
                        #hoverlabel=dict(bgcolor='rgba(255, 255, 255, 0.8)',
                        #                bordercolor='black',
                        #                align='left',
                        #                namelength=-1),
                        #hovertext=hover_text,
                        showlegend=showlegend,
                    ))
                    
                    fig.add_annotation(
                        x=x_position + text_offset,
                        y=y_pos,
                        text=FullName2AbbrName.get(feature.name, feature.name),
                        showarrow=False,
                        font=dict(size=gene_label_size, color=gene_label_color),
                        xanchor='left',
                        yanchor='middle', 
                        hoverlabel=dict(bgcolor='rgba(255, 255, 255, 0.8)',
                                        bordercolor='black'),
                        hovertext=hover_text,
                    )
                    
                    if i+1 < len(features[1:]):
                        if features[1:][i+1].location.strand == -1:                            
                            x_position += self._get_box_param(FullName2AbbrName.get(features[1:][i+1].name, features[1:][i+1].name))[2] + self._get_box_param( FullName2AbbrName.get(features[1:][i].name, features[1:][i].name))[2]
                        else:
                            x_position += self._get_box_param(FullName2AbbrName.get(features[1:][i].name, features[1:][i].name))[2]
                
                elif feature.location.strand == -1:
                    fig.add_trace(go.Scatter(
                        x=[x_position, x_position - tail_length, x_position - tail_length - head_length,
                          x_position - tail_length, x_position, x_position],
                        y=[y_pos - height/2,
                           y_pos - height/2,
                           y_pos,
                           y_pos + height/2,
                           y_pos + height/2,
                           y_pos - height/2],
                        mode='lines',
                        fill="toself",
                        fillcolor=fill_color,
                        line=dict(color='black', width=1),
                        name=group_name,
                        legendgroup=group_name,
                        hoverinfo='skip',
                        #hoverinfo='text',
                        #hoveron='points+fills',
                        #hoverlabel=dict(bgcolor='rgba(255, 255, 255, 0.8)',
                        #                bordercolor='black',
                        #                align='left',
                        #                namelength=-1),
                        #hovertext=hover_text,
                        showlegend=showlegend,
                    ))
                    
                    fig.add_annotation(
                        x=x_position - text_offset,
                        y=y_pos,
                        text=FullName2AbbrName.get(feature.name, feature.name),
                        showarrow=False,
                        font=dict(size=gene_label_size, color=gene_label_color),
                        xanchor='right',
                        yanchor='middle',
                        hoverlabel=dict(bgcolor='rgba(255, 255, 255, 0.8)',
                                        bordercolor='black'),
                        hovertext=hover_text,
                    )
                    
                    if i+1 < len(features[1:]):
                        if features[1:][i+1].location.strand == -1:
                            x_position +=  self._get_box_param(FullName2AbbrName.get(features[1:][i+1].name, features[1:][i+1].name))[2]
                            
                elif feature.location.strand == 0:
                    _staus_brake = True
                    fig.add_trace(go.Scatter(
                        x=[x_position+box_width*1/6, x_position+box_width/2, ],
                        y=[y_pos-height/2, y_pos + height/2, ],
                        mode='lines',
                        fill="toself",
                        line=dict(color='black', width=1),
                        showlegend=False,
                    ))
                    fig.add_trace(go.Scatter(
                        x=[x_position+box_width/2, x_position+box_width-box_width*1/6],
                        y=[y_pos-height/2, y_pos + height/2, ],
                        mode='lines',
                        fill="toself",
                        line=dict(color='black', width=1),
                        showlegend=False,
                    ))
                    #fig.add_trace(go.Scatter(
                    #    x=[x_position, x_position+box_width],
                    #    y=[y_pos, y_pos],
                    #    mode='lines',
                    #    fill="toself",
                    #    line=dict(color='black', width=1),
                    #    showlegend=False,
                    #))
                    
                    if i+1 < len(features[1:]):
                        if features[1:][i+1].location.strand == -1:
                            x_position +=  self._get_box_param(FullName2AbbrName.get(features[1:][i+1].name, features[1:][i+1].name))[2] + self._get_box_param(FullName2AbbrName.get(features[1:][i].name, features[1:][i].name))[2]
                        else:
                            x_position += self._get_box_param(FullName2AbbrName.get(features[1:][i].name, features[1:][i].name))[2]

                            
                tmp_groups.add(group_name)
                
        layout_dict = {#'template':'ggplot2',
            'font':  dict(size=12),
            'xaxis': dict(showticklabels=False, showgrid=False, zeroline=False),
            'yaxis': dict(tickvals=y_positions,
                          ticktext=species_labels,
                          fixedrange=False,
                          constrain='domain',
                          tickfont=dict(size=species_label_size, style='italic', color=species_label_color)
            ),
            'showlegend': show_legend,
            'margin':dict(l=50, r=50, t=50, b=50, pad=10),
            'height': 175 + 50*len(genomes),
            'width': 1500,
            'legend': dict(x=1, y=0, xanchor="right", yanchor="top", traceorder='reversed',
                           bgcolor="rgba(255, 255, 255, 0.5)",orientation="h",font=dict(color='black'),
                           bordercolor="black", borderwidth=1)
        }

        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(symbol='triangle-left-open', color='black', size=12),
            name='Reverse',
            showlegend=True
        ))

        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(symbol='triangle-right-open', color='black', size=12),
            name='Forward',
            showlegend=True
        ))
        
        if  _staus_brake:
            fig.add_trace(go.Scatter(
                x=[None], y=[None],
                mode='markers',
                marker=dict(symbol='line-ne', color='black'),
                name='// Break',
                showlegend=True
            ))
            
        layout_dict.update(self.default_layout)
        fig.update_yaxes(range=[-1, len(genomes)+1])
        
        fig.update_layout(**layout_dict) 
        if output:
            fig.write_html(output, config=self.nav_config)
        #else:
        return fig
            
    def draw_linear_MT_plotly(self,
                              files,
                              output=None,
                              abbr=False,
                              isfilename2species=False,
                              colors="mitofish",
                              show_gene_label=True,
                              gene_label_size=9,
                              gene_label_color='black',
                              species_label_size=16,
                              species_label_color='black',
                              show_legend=True,
                              start=None,
                              tidyname=False,
                              add_id=False,
                              editable=False,
                              force_reoriented=False,
                              remove_NCR = False,
                             ):
        """
        Descripton:
            The order of mitochondrial genes, arranged in proportion to their size.
            Interactive chart developed using Plotly.
        
        Parameters:
            files: {str, list, tuple} one or more genbankfile or NCBI accession ID.
            output: {str} a path of fig save, save as HTML file.
            abbr: {bool} whether to abbreviate species names.
            isfilename2species: {bool} whether filename convert to species.
            colors: {str, dict} themes such as, Chen, Tan, ogdraw, mitofish,
                                mitofish1, mitoz,  gggenes, chloroplot, grey, igv.
            gene_label_size: {int} gene label size.
            gene_label_color: {str} gene label color.
            species_label_size: {int} species name size.
            species_label_color: {int} species name color.
            show_legend: {bool} show fig legend.
            start: {None, str} initial feature, such as, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
                         COX1, COX2, COX3, ATPase6, ATPase8, Cytb, tRNA-His, tRNA-Pro,
                         tRNA-Thr, tRNA-Trp, tRNA-Met, tRNA-Asp, tRNA-Ala, tRNA-Gln,
                         tRNA-Ile, tRNA-Arg, tRNA-Tyr, tRNA-Phe, tRNA-Lys, tRNA-Gly,
                         tRNA-Asn, tRNA-Leu, tRNA-Glu, tRNA-Val, tRNA-Cys, tRNA-Ser,
                         12S rRNA, 16S rRNA, D-loop. default=None.
            add_id: {bool} Species add to accession id from NCBI.
            editable: {bool} Make HTML files editable.
            force_reoriented: {bool} force-reoriendted linear mtgenome.
            remove_NCR {bool}: Do not display D-loop; this may be useful for comparing a mix of genbank with and without NCR annotation.     
        """
        
        color2groups = self._get_colors_group(colors)
        tmp_groups = set()
        
        if not isinstance(files, (list, tuple)):
            files = [files]

        if editable:
            self.nav_config["editable"] = True
            
        subplot_height = 100
        fixed_spacing = 50
        total_height = (subplot_height * (len(files)+1)) + (fixed_spacing * ((len(files)+1) - 1))
        
        fig = make_subplots(rows=len(files)+1, cols=1, subplot_titles=[], vertical_spacing=fixed_spacing / total_height,
                            shared_xaxes=True)
        
        genome_lengths = []
        species_names = []
        _staus_brake = False
        for i, file in enumerate(files):
            features = get_features(file, abbr=abbr, isfilename2species=isfilename2species, colors=colors, start=start, force_reoriented=force_reoriented)
            
            if remove_NCR:
                features = [feature for feature in features if feature.name != "D-loop"]
            
            #if tidyname:
            #    for j, f in enumerate(features):
            #        features[j].name = FullName2AbbrName.get(f.name, f.name)
            
            genome_length = len(features[0].location)
            genome_lengths.append(genome_length)
            row = i + 1
            
            medial_axis = 0
            species_name = features[0].name + (f" ({features[0].accession})" if add_id else "")
            species_names.append(species_name)
            
            fig.add_trace(go.Scatter(
                x=[0, genome_length, genome_length, 0, 0],
                y=[medial_axis, medial_axis, medial_axis, medial_axis, medial_axis],
                fill="toself",
                mode='lines',
                name=features[0].name,
                fillcolor='black',
                line=dict(color='black', width=1),
                showlegend=False
            ), row=row, col=1)
            _staus_brake_tmp = False
            for index, feature in enumerate(features[1:]):    
                start_pos = int(feature.location.start)
                end_pos = int(feature.location.end)                
                #length = end_pos - start_pos + 1
                fill_color = self._convert_color(feature.color)
                hover_text = f"""Gene Name: {feature.name}<br>Position: {start_pos+1}-{end_pos}{'(+) ' if feature.location.strand == 1 else '(-) '}"""
                if _staus_brake_tmp:
                    start_pos += 200
                    end_pos += 200
                
                #if color2groups[feature.name].get('group', "Other") not in tmp_groups:
                #    showlegend = True
                #else:
                #    showlegend = False
                    
                group_name = color2groups[feature.name]['group'] if feature.name in color2groups else color2groups["Other genes"]['group']
                if group_name not in tmp_groups:
                    showlegend = True
                else:
                    showlegend = False
                    
                if feature.location.strand == 1:
                    if show_gene_label:
                        fig.add_annotation(
                            x=(start_pos + end_pos) / 2,
                            y=medial_axis+0.55,
                            #text=feature.name,
                            text=FullName2AbbrName.get(feature.name, feature.name) if tidyname else feature.name,
                            
                            showarrow=False,
                            textangle=90,
                            font=dict(size=gene_label_size, color=gene_label_color),
                            xanchor='center',
                            yanchor='bottom',
                            hovertext=hover_text,
                            hoverlabel=dict(bgcolor='rgba(255, 255, 255, 0.8)', bordercolor='black', 
                                            ),
                            row=row,
                            col=1,
                        )
                        
                    fig.add_trace(go.Scatter(
                        x=[start_pos, end_pos, end_pos, start_pos, start_pos],
                        y=[medial_axis+0.5, medial_axis+0.5, medial_axis, medial_axis, medial_axis+0.5],
                        fill="toself",
                        mode='lines',
                        legendgroup=group_name,
                        name=group_name,
                        fillcolor=fill_color,
                        line=dict(color='black', width=1),
                        hoverinfo='skip',
                        #hoverinfo='text',
                        #hoveron='points+fills',
                        #hoverlabel=dict(bgcolor='rgba(255, 255, 255, 0.8)',
                        #                bordercolor='black',
                        #                align='left',
                        #                namelength=-1),
                        #hovertext=hover_text,
                        showlegend=showlegend), row=row, col=1)    
                    
                    fig.add_trace(go.Scatter(
                        x=[(end_pos - start_pos)/2 + start_pos],
                        y=[medial_axis+0.25],
                        #fill="toself",
                        mode='lines',
                        line=dict(color=None, width=1),
                        hoverinfo='text',
                        hoveron='points+fills',
                        hoverlabel=dict(bgcolor='rgba(255, 255, 255, 0.8)',
                                        bordercolor='black',
                                        align='left',
                                        namelength=-1),
                        hovertext=hover_text,
                        showlegend=False
                    ), row=row, col=1)
                    
                elif feature.location.strand == -1:
                    if show_gene_label:
                        
                        fig.add_annotation(
                            x=(start_pos + end_pos) / 2,
                            y=medial_axis-0.55,
                            #text=feature.name,
                            text=FullName2AbbrName.get(feature.name, feature.name) if tidyname else feature.name,
                            showarrow=False,
                            textangle=90,
                            font=dict(size=gene_label_size, color=gene_label_color),
                            hovertext=hover_text,
                            hoverlabel=dict(bgcolor='rgba(255, 255, 255, 0.8)', bordercolor='black'),
                            xanchor='center',
                            yanchor='top',
                            row=row,
                            col=1
                        )
                        
                    fig.add_trace(go.Scatter(
                        x=[start_pos, end_pos, end_pos, start_pos, start_pos],
                        y=[medial_axis, medial_axis, medial_axis-0.5, medial_axis-0.5, medial_axis],
                        fill="toself",
                        mode='lines',
                        name=group_name,
                        fillcolor=fill_color,
                        line=dict(color='black', width=1),
                        legendgroup=group_name,
                        hoverinfo='skip',
                        #hoverinfo='text',
                        #hoveron='points+fills',
                        #hoverlabel=dict(bgcolor='rgba(255, 255, 255, 0.8)',
                        #                bordercolor='black',
                        #                align='left',
                        #                namelength=-1),
                        #hovertext=hover_text,
                        showlegend=showlegend
                    ), row=row, col=1)
                    
                    fig.add_trace(go.Scatter(
                        x=[(end_pos - start_pos)/2 + start_pos],
                        y=[medial_axis-0.25],
                        #fill="toself",
                        mode='lines',
                        line=dict(color=None, width=1),
                        hoverinfo='text',
                        hoveron='points+fills',
                        hoverlabel=dict(bgcolor='rgba(255, 255, 255, 0.8)',
                                        bordercolor='black',
                                        align='left',
                                        namelength=-1),
                        hovertext=hover_text,
                        showlegend=False
                    ), row=row, col=1)

                else:
                    _staus_brake = True
                    _staus_brake_tmp = True
                    box_width = 200
                    x_position = end_pos
                    #x_position = int(features[index].location.end)
                    
                    fig.add_trace(go.Scatter(
                        x=[genome_length, genome_length+box_width, genome_length+box_width, genome_length, genome_length],
                        y=[medial_axis, medial_axis, medial_axis, medial_axis, medial_axis],
                        fill="toself",
                        mode='lines',
                        name=features[0].name,
                        fillcolor='black',
                        line=dict(color='black', width=1),
                        showlegend=False
                    ), row=row, col=1)
                    
                    fig.add_trace(go.Scatter(
                        x=[x_position+4, x_position+box_width],
                        y=[medial_axis, medial_axis],
                        mode='lines',
                        fill="toself",
                        line=dict(color='white', width=3),
                        showlegend=False,
                    ), row=row, col=1)
                    
                    fig.add_trace(go.Scatter(
                        x=[x_position+box_width*1/6, x_position+box_width/2],
                        y=[medial_axis-0.5, medial_axis + 0.5],
                        mode='lines',
                        fill="toself",
                        line=dict(color='black', width=1),
                        showlegend=False,
                    ), row=row, col=1)
                    fig.add_trace(go.Scatter(
                        x=[x_position+box_width/2, x_position+box_width-box_width*1/6],
                        y=[medial_axis-0.5, medial_axis + 0.5],
                        mode='lines',
                        fill="toself",
                        line=dict(color='black', width=1),
                        showlegend=False,
                    ), row=row, col=1)
		    
                tmp_groups.add(group_name)
                
        for i in range(1, len(files)+1):
            fig.update_yaxes(range=[-1, 1],
                             tickvals=[medial_axis],ticktext=[species_names[i-1]],
                             tickfont=dict(size=species_label_size, style='italic', color=species_label_color),
                             row=i, col=1, showticklabels=True)
            if i < len(files):
                fig.update_xaxes(
                    showline=False,
                    showticklabels=False,
                    showgrid=False,
                    row=i, col=1
                )
            else:
                fig.update_xaxes(range=[0, max(genome_lengths)+500],
                    linewidth=1,
                    showline=True,
                    linecolor="black",
                    showticklabels=True,
                    showgrid=False,
                    row=i, col=1,
                    tickfont=dict(color='black'),
                )

        if  _staus_brake:
            fig.add_trace(go.Scatter(
                x=[None], y=[None],
                mode='markers',
                marker=dict(symbol='line-ne', color='black'),
                name='// Break',
                showlegend=True
            ))
            
        layout_dict = {'showlegend': show_legend,
                       'width': 2000,
                       'height': 100+total_height,
                       'margin': dict(l=50, r=50, t=50, b=50, pad=40),
                       'legend': dict(x=1, y=0, xanchor="right", yanchor="bottom", traceorder='reversed',
                                      bgcolor="rgba(255, 255, 255, 0.5)",orientation="h",font=dict(color='black'),
                                      bordercolor="black", borderwidth=1)
                      }
        layout_dict.update(self.default_layout)
        fig.update_layout(**layout_dict)
        
        if output:
            fig.write_html(output, config=self.nav_config)
        #else:
        return fig
            
            
def draw_linear_MT_nonproportional_interactive(files, output=None, abbr=False, isfilename2species=False,
                                               colors="mitofish", gene_label_size=9, gene_label_color='black',
                                               show_legend=True,
                                               start=None, species_label_size=12, species_label_color='black', add_id=False, editable=False, force_reoriented=False, remove_NCR = False):
    """
    Descripton:
        The order of mitochondrial genes is not depicted in proportion to their gene size.
        Interactive chart developed using Plotly.
    
    Parameters：
        files: {str, list, tuple} one or more genbankfile or NCBI accession ID.
        output: {str} a path of fig save, save as HTML file.
        abbr: {bool} whether to abbreviate species names.
        isfilename2species: {bool} whether filename convert to species.
        colors: {str, dict} themes such as, Chen, Tan, ogdraw, mitofish,
                            mitofish1, mitoz,  gggenes, chloroplot, grey, igv.
        gene_label_size: {int} gene label size.
        gene_label_color: {str} gene label color.
        species_label_size: {int} species name size.
        species_label_color: {str} species name color.
        show_legend: {bool} show fig legend.
        start: {None, str} initial feature, such as, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
                     COX1, COX2, COX3, ATPase6, ATPase8, Cytb, tRNA-His, tRNA-Pro,
                     tRNA-Thr, tRNA-Trp, tRNA-Met, tRNA-Asp, tRNA-Ala, tRNA-Gln,
                     tRNA-Ile, tRNA-Arg, tRNA-Tyr, tRNA-Phe, tRNA-Lys, tRNA-Gly,
                     tRNA-Asn, tRNA-Leu, tRNA-Glu, tRNA-Val, tRNA-Cys, tRNA-Ser,
                     12S rRNA, 16S rRNA, D-loop. default=None.
        add_id: {bool} Species add to accession id from NCBI.
        editable: {bool} Make HTML files editable.
        force_reoriented: {bool} force-reoriendted linear mtgenome.
        remove_NCR {bool}: Do not display D-loop; this may be useful for comparing a mix of genbank with and without NCR annotation. 
    """
    return InteractiveMTVisualizer().draw_linear_MT_nonproportional_plotly(files=files,
                                                                           output=output,
                                                                           abbr=abbr,
                                                                           isfilename2species=isfilename2species, 
                                                                           colors=colors,
                                                                           gene_label_size=gene_label_size, 
                                                                           gene_label_color='black',
                                                                           show_legend=show_legend,
                                                                           start=start,
                                                                           species_label_size=species_label_size,
                                                                           species_label_color=species_label_color,
                                                                           add_id=add_id,
                                                                           editable=editable,
                                                                           force_reoriented=force_reoriented,
                                                                           remove_NCR =remove_NCR)

def draw_linear_MT_interactive(files, output=None,
                               abbr=False,
                               isfilename2species=False,
                               colors="mitofish",
                               show_gene_label=True,
                               gene_label_size=9,
                               gene_label_color='black',
                               species_label_size=16,
                               species_label_color='black',
                               show_legend=True,
                               start=None,
                               tidyname=False,
                               add_id=False,
                               editable=False,
                               force_reoriented=False,
                               remove_NCR = False):
    """
    Descripton:
        The order of mitochondrial genes, arranged in proportion to their size.
        Interactive chart developed using Plotly.
    
    Parameters：
        files: {str, list, tuple} one or more genbankfile or NCBI accession ID.
        output: {str} a path of fig save, save as HTML file.
        abbr: {bool} whether to abbreviate species names.
        isfilename2species: {bool} whether filename convert to species.
        colors: {str, dict} themes such as, Chen, Tan, ogdraw, mitofish,
                            mitofish1, mitoz,  gggenes, chloroplot, grey, igv.
        gene_label_size: {int} gene label size.
        gene_label_color: {str} gene label color.
        species_label_size: {int} species name size.
        species_label_color: {str} species name color.
        show_legend: {bool} show fig legend.
        start: {None, str} initial feature, such as, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
                     COX1, COX2, COX3, ATPase6, ATPase8, Cytb, tRNA-His, tRNA-Pro,
                     tRNA-Thr, tRNA-Trp, tRNA-Met, tRNA-Asp, tRNA-Ala, tRNA-Gln,
                     tRNA-Ile, tRNA-Arg, tRNA-Tyr, tRNA-Phe, tRNA-Lys, tRNA-Gly,
                     tRNA-Asn, tRNA-Leu, tRNA-Glu, tRNA-Val, tRNA-Cys, tRNA-Ser,
                     12S rRNA, 16S rRNA, D-loop. default=None.
        add_id: {bool} Species add to accession id from NCBI.
        editable: {bool} Make HTML files editable.
        force_reoriented: {bool} force-reoriendted linear mtgenome.
        remove_NCR {bool}: Do not display D-loop; this may be useful for comparing a mix of genbank with and without NCR annotation.
    """
    return InteractiveMTVisualizer().draw_linear_MT_plotly(files=files,
                                                           output=output,
                                                           abbr=abbr,
                                                           isfilename2species=isfilename2species,
                                                           colors=colors,
                                                           show_gene_label=show_gene_label,
                                                           gene_label_size=gene_label_size,
                                                           gene_label_color=gene_label_color,
                                                           species_label_size=species_label_size,
                                                           species_label_color=species_label_color,
                                                           show_legend=show_legend,
                                                           start=start,
                                                           tidyname=tidyname,
                                                           add_id=add_id,
                                                           editable=editable,
                                                           force_reoriented=force_reoriented,
                                                           remove_NCR =remove_NCR)
