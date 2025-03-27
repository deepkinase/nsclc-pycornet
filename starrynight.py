# -*- coding: utf-8 -*-

"""
Script Name: StarryNight
Description: Illuminating the Kinase Activity Starmap in NSCLC using tyrosine phosphoproteomics data
Author: HuiZhou
Date: 2024-07-28
"""

import pandas as pd
import json,os,sys,re
import statistics,subprocess
import math,re,base64,argparse
desc = "starrynight for bx project 20240130"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-prefix",help="prefix")
parser.add_argument("-files",help="files of inputs")
parser.add_argument("-cates",help="category for features")
parser.add_argument("-neg",help="cutoff for negative rho")
parser.add_argument("-pos",help="cutoff for positive rho")
parser.add_argument("-pval",help="cutoff for pvalue")
parser.add_argument("-outdir",help="outdir")
args = parser.parse_args()

prefix = args.prefix
files = args.files
cates = args.cates
neg = args.neg
pos = args.pos
pval = args.pval
outdir = args.outdir
names = outdir+'/'+prefix 
###### read data
net_all = pd.read_csv(files)
cates_data =  pd.read_csv(cates)
pro_data = cates_data[cates_data['Source']=='Proteins']['name'].tolist()
pros = pro_data
network = {}
node_list = []
edge_list = []
###### Screening the overall data based on the parameters passed.
net_all = net_all[((net_all['Spearman_Rho'] >= float(pos)) | (net_all['Spearman_Rho'] <= float(neg))) & (net_all['Spearman_Pvalue'] <= float(pval))]
net_cc =  net_all[net_all['Category']=='Compound.vs.Compound'].sort_values(by='Spearman_Pvalue', ascending=True)
net_gg = net_all[net_all['Category']=='Protein.vs.Protein'].sort_values(by='Spearman_Pvalue', ascending=True)
###### data_url
symbol_path = './fig.svg.motify/'
def url_fun(svg):
    with open(svg, "r") as f:
        svg_content = f.read()
    base64_data = base64.b64encode(svg_content.encode()).decode()
    data_url = f"data:image/svg+xml;base64,{base64_data}"
    return(data_url)
star4p = url_fun(symbol_path+'star4p.svg')
star4 = url_fun(symbol_path+'star4.svg')
benzene = url_fun(symbol_path+'benzene.svg')
url = {'benzene':benzene,'star4':star4,'star4p':star4p}

###### Replace the content of the legend.
def find_legend(file_path, keyword):
    files = open(file_path,'r').readlines()
    for i, line in enumerate(files):
        if keyword in line:
            return(i+3)
    return matching_lines

###### Define the node size based on the rank of the p-value.
###### The weight of degree (connectivity) on node size.
def degree(df):
    return(df.shape[0])

def text_size(num):
    if num >= 100 :
        return(80)
    else:
        return(num)

def p_node(net_gc,fe):
    net_gc.loc[:,'rank'] = net_gc['Spearman_Pvalue'].rank(ascending=False)
    net_gc['rank_cubed'] = net_gc['rank'].apply(lambda x: x ** 3)
    net_gc = net_gc.sort_values(by='rank_cubed',ascending=False)
    net_gc.to_csv(names+'net_rank.csv')
    fe_gc = net_gc['Feature1'].tolist()+net_gc['Feature2'].tolist()
    fe_sele = net_gc[(net_gc['Feature1'] == fe) | (net_gc['Feature2'] == fe)]
    p_att = abs(math.log10(fe_sele.iloc[0]['Spearman_Pvalue']+1e-10))
    rho_att = abs(fe_sele.iloc[0]['Spearman_Rho'])
    de_att = degree(fe_sele)
    size_node = round(de_att*1.5+4,2)

    return(size_node)

###### Define the node size.
def node_size(net_gc,fe):
    fe_gc = net_gc['Feature1'].tolist()+net_gc['Feature2'].tolist()
    if fe not in fe_gc:
        nodesize = 4
    else:
        nodesize = p_node(net_gc,fe)
    return(nodesize)
###### change name 
def new_name(name):
    pattern = r'\((.*?)\)'
    if '(' in name:
        matches = re.findall(pattern, name)[0]
        matches = name.split('(')[0].strip()+'_'+matches

    else:
        matches = name.strip()
    return(matches)

def names_fun(name):
    if '(' in name:
        return(name.split('(')[0].strip())

    elif len(re.findall(r'(.+)_\d+',name)) == 1:
        return(re.findall(r'(.+)_\d+',name)[0])
#        return

    else:
        return(name.strip())
###### change drug id_name to drug name
def network_extend(net_gc,number):
    for n in list(range(number)):
        if n != number:
            features_net_gc = net_gc['Feature1'].tolist()+net_gc['Feature2'].tolist()
            net_cc_sele_extend = net_cc[~(net_cc['Feature1'].isin(features_net_gc) | net_cc['Feature2'].isin(features_net_gc))]
            net_gg_sele_extend = net_gg[~(net_gg['Feature1'].isin(features_net_gc) | net_gg['Feature2'].isin(features_net_gc))]
            net_cc_gg_sele_extend = pd.concat([net_cc_sele_extend,net_gg_sele_extend])
        net_gc = pd.concat([net_gc,net_cc_gg_sele_extend])
    net_gc = net_gc.drop_duplicates()
    return(net_gc)

net_gc_extend =  net_gg

###### The Gene-Compounds network serves as the core network, extended by the GG and CC networks.
net_sele = net_gc_extend
###### Filter out links with low degree, retaining features with a degree of at least 2.
features = net_sele['Feature1'].tolist()+net_sele['Feature2'].tolist()
net_list  = [net_sele.iloc[n].tolist() for n in list(range(net_sele.shape[0])) if features.count(net_sele.iloc[n]['Feature1']) >=2 and features.count(net_sele.iloc[n]['Feature2']) >=2]
net = pd.DataFrame(net_list, columns=net_sele.columns)

def net_filt_fun(net):
    net_gg = net[net['Category'] == 'Protein.vs.Protein']
    net_gg_n = net_gg.shape[0]

    if net_gg.shape[0] < 1000:
        net_final = net_gg

    else:
        net_gg = net_gg.sort_values(by='Spearman_Pvalue').head(1000)
        net_final = net_gg
    return(net_final)

net = net_filt_fun(net)
###### Output the final network.
net.to_csv(names+'.GC.extend.add.gg.cc.'+'csv')
###### Features category
features_list = net['Feature1'].tolist()+net['Feature2'].tolist()
features_set_change = [ names_fun(fe) for fe in set(features_list)]

###### node_size_distribution
node_size_distr = []
for fe in set(features_list):
    size = features_list.count(fe)
    cate = cates_data[cates_data['name']==fe]['category'].tolist()
    name = names_fun(fe) if features_set_change.count(names_fun(fe)) == 1 else new_name(fe)
    shape = cates_data[cates_data['name']==fe]['Symbol'].tolist()
    colour = cates_data[cates_data['name']==fe]['Color'].tolist()
    symbolsize = node_size(net,fe)
    node_size_distr.append(symbolsize)
    node_list.append({"name":name,"symbolSize": symbolsize,"value": symbolsize, "category": cate[0], "symbol":'image://'+url[shape[0]],"draggable" : "True","selected":"True","focusnode":"True","roam":"True"})
###### links
for n  in list(range(net.shape[0])):
    node1 = net.iloc[n]['Feature1']
    node2 = net.iloc[n]['Feature2']
    node1 = names_fun(node1) if features_set_change.count(names_fun(node1)) == 1 else new_name(node1)
    node2 = names_fun(node2) if features_set_change.count(names_fun(node2)) == 1 else new_name(node2)
    value = net.iloc[n]['link_cate'] #20240311
    edge_list.append({"source":node1,"target":node2,"value":value})
###### category
categories_lists = ['Biomarker']
categories = [{"name":ca} for ca in categories_lists]
###### network json
network = {
    "nodes": node_list,
    "links": edge_list,
}
###### Write JSON data to a file
output_filename = names+'.jason'
with open(output_filename, "w") as outfile:
    json.dump(network, outfile, indent=4)

###### Draw the graph
from pyecharts import options as opts
from pyecharts.charts import Graph
nodes = node_list
links = edge_list
outhtml = names+'.GC.extend.add.gg.cc.'+'html'
###### Display the top 10 symbol sizes at most.
node_size_distr_sort = sorted(node_size_distr, reverse=True)
for node in nodes:
    if node["symbolSize"] > node_size_distr_sort[10]:
        if node["symbolSize"] > 120:
            node["label"] = {"show": False,"fontSize":30}
        else:
            node["label"] = {"show": False,"fontSize": node["value"]/2.5}
    else:
        node["label"] = {"show": False}


color_mapping = {
    0:"#E0E0E0",
    7:"#E0E0E0",
    21:"#E0E0E0",
    35:"#E0E0E0",
    1 : "#FF6666",
    9 : "#3399FF",
    25 : "#009900",
    49 : "#FF8000",
    3 : "#CC00CC",
    5 : "#999900",
    15 :"#00CCCC",
}

lines = []
for link in links:
    line = {
        "source": link["source"],
        "target": link["target"],
        "value": link["value"],
        "lineStyle": {
            "color":  color_mapping[link["value"]],
            "width": 2.5,
            "type":  "dashed",
            "opacity": 0.5,
        },
    }
    lines.append(line)
###### Add an external stylesheet link.
link_html = """
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Exo:wght@500&display=swap" rel="stylesheet">
"""
###### Customize colors.
custom_colors = [
    "#fbefcb",
    "#91cc75",
    "#fac858",
    "#ee6666",
    "#73c0de",
    "#3ba272",
    "#fc8452",
    "#9a60b4",
    "#ea7ccc",
]
graph = (
        Graph(init_opts=opts.InitOpts(
            width="100vw", 
            height="100vh",
            page_title="DeepKinase StarryNight",
            theme = "white",
            animation_opts=opts.AnimationOpts(animation=True,
                 animation_duration=40000,
                 animation_easing="cubicInOut",
                 animation_delay=1000,
                 animation_duration_update=40000,
                 animation_easing_update="cubicInOut",
                 animation_delay_update=1000,
                 animation_threshold=30000,

                 ),
                )
            )
        .add("", 
            nodes,
            lines,
            categories=categories,
            repulsion=800,
            gravity=0.2,
            friction=0.1,
            is_rotate_label=True,
            linestyle_opts=opts.LineStyleOpts(curve=0.2),
            layout="force",
            label_opts=opts.LabelOpts(position="right",font_family="Exo"),
        )
        .set_global_opts(
            title_opts=opts.TitleOpts(title="DK StarryNight: "+names.split('/')[-1].replace('_',' '),
                padding = 30,
                pos_top="top", 
                pos_left="center",
                title_textstyle_opts=opts.TextStyleOpts(
                     color= "#C0C0C0",
                     font_family= "Exo",
                     font_size = 30,
                     )
                ),
            legend_opts=opts.LegendOpts(selected_mode="multiple",orient="vertical",pos_left="2%", pos_top="3%",
                    textstyle_opts=opts.TextStyleOpts(
                     color= "#fbefcb99",
                     font_family= "Exo",
                     font_size = 12,
                     )
),  # Allow selecting multiple categories

            )
#        .set_colors(custom_colors)
    .render(outhtml)
)
###### Add Google Fonts and replace the legend.
f1 = './drug.cate_pro.cate/google.style.html'
f2 = './drug.cate_pro.cate/legends.html'
f3 = './drug.cate_pro.cate/legend.html'
f4 = outhtml
legend_index = find_legend(f4,"tooltip")
command1 = f"sed -i '7r {f1}' {f4}"
command2 = f'sed -i "s/<body >/<body style=\\"margin:0\\">/g" {f4}'
command3 = f"sed -i '{legend_index}r {f2}' {f4}"
command4 = f"sed -i '/{f3}/d' {f4}"
os.system(command1)
os.system(command2)
os.system(command3)
