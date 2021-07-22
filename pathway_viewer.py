#!/usr/bin/env python

from PyQt5.QtWidgets import *
from PyQt5.QtSvg import *
from PyQt5.QtCore import *
from PyQt5.QtWebEngineWidgets import *
from network_builder import *
import pandas as pd
import networkx as nx

def make_coexpression_matrix(expression_db, gene_db):
    gpl570 = pd.read_csv(gene_db, sep="\t", low_memory=False)[
        ["ID", "Gene Symbol"]
    ]
    gse = pd.read_csv(expression_db, sep="\t", low_memory=False)

    tribe2_seq = pd.read_csv("TRIBE2_seq_res.csv")
    dataset_genes = tribe2_seq["Biomarker"].unique()

    gpl570 = gpl570[gpl570["Gene Symbol"].isin(dataset_genes)].set_index("ID")
    gene_data = gpl570.join(gse.set_index("ID_REF"))
    gene_expression = gene_data.groupby("Gene Symbol").mean()

    coexpression = gene_expression.T.corr()
    return coexpression

coex = make_coexpression_matrix("GSE40367-stripped.txt", "GPL570-stripped.txt")
active_pathway = None

@pyqtSlot()
def export_gexf():
    path = QFileDialog.getSaveFileName()
    nx.write_gexf(active_pathway.graph, path)

app = QApplication([])
svg_pane = QWebEngineView()
def svg_changed(svg_data):
    print('setting content')
    svg_pane.setContent(svg_data, "image/svg+xml")
    print('content set')

def make_graph(threshold):
    active_pathway = make_pathway_from_thres(threshold, coex)
    print('Made pathway', threshold)
    dot_graph = nx.nx_pydot.to_pydot(active_pathway.graph)
    print('Made graph')
    svg = dot_graph.create_svg()
    print('Made svg')
    svg_changed(svg)


window = QWidget()

right_pane = QWidget()
thrs_box = QDoubleSpinBox()
thrs_box.setMinimum(0)
thrs_box.setMaximum(1.0)
thrs_box.setSingleStep(0.05)
thrs_box.valueChanged.connect(make_graph)
thrs_label = QLabel('Link threshold')
export_button = QPushButton('Save as GEXF')
export_button.clicked.connect(export_gexf)

right_layout = QVBoxLayout()
right_layout.addWidget(thrs_label)
right_layout.addWidget(thrs_box)
right_layout.addWidget(export_button)
right_pane.setLayout(right_layout)

mainlayout = QGridLayout()
mainlayout.addWidget(svg_pane, 0, 0, 1, 2)
mainlayout.addWidget(right_pane, 0, 2, 1, 1)

window.setLayout(mainlayout)
window.show()

app.exec()
