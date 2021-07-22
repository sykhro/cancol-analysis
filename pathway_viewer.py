#!/usr/bin/env python

from PyQt5.QtWidgets import *
from PyQt5.QtSvg import *
from PyQt5.QtCore import *
from PyQt5.QtWebEngineWidgets import *
from network_builder import *
import pandas as pd
import networkx as nx
import sys

DEFAULT_EXP = "GSE40367-stripped.txt"
DEFAULT_GENES = "GPL570-stripped.txt"

def make_coexpression_matrix(expression_db, gene_db):
    gpl570 = pd.read_csv(gene_db, sep="\t", low_memory=False)[["ID", "Gene Symbol"]]
    gse = pd.read_csv(expression_db, sep="\t", low_memory=False)

    tribe2_seq = pd.read_csv("TRIBE2_seq_res.csv")
    dataset_genes = tribe2_seq["Biomarker"].unique()

    gpl570 = gpl570[gpl570["Gene Symbol"].isin(dataset_genes)].set_index("ID")
    gene_data = gpl570.join(gse.set_index("ID_REF"))
    gene_expression = gene_data.groupby("Gene Symbol").mean()

    coexpression = gene_expression.T.corr()
    return coexpression

class PathwayView(QWidget):
    def __init__(self, pathway):
        super().__init__()
        self.active_pathway = pathway

class Viewer(QWidget):
    def __init__(self):
        super().__init__()

        self.active_pathway = None
        self.exp_path = DEFAULT_EXP
        self.genes_path = DEFAULT_GENES
        self.ref_coex = make_coexpression_matrix(self.exp_path, self.genes_path)

        # Settings -----------------------------------------------------------------------------
        settings_group = QGroupBox("Pathway settings")
        thrs_layout = QHBoxLayout()
        thrs_box = QDoubleSpinBox()
        thrs_box.setMinimum(0)
        thrs_box.setMaximum(1.0)
        thrs_box.setSingleStep(0.05)
        thrs_box.valueChanged.connect(self.make_graph)
        thrs_label = QLabel("Link threshold")
        thrs_layout.addWidget(thrs_label)
        thrs_layout.addWidget(thrs_box)

        genes_path = QLineEdit()
        genes_path.setReadOnly(True)
        genes_path.setText(self.genes_path)
        genes_selection = QPushButton("Genes")
        genes_selection.clicked.connect(self.set_genes_db)

        expression_path = QLineEdit()
        expression_path.setReadOnly(True)
        expression_path.setText(self.exp_path)
        expression_selection = QPushButton("Expression data")
        expression_selection.clicked.connect(self.set_expression_data)
        
        settings_layout = QGridLayout()
        settings_layout.addWidget(thrs_label, 0, 0)
        settings_layout.addWidget(thrs_box, 0, 1)
        warn = QLabel("Warning! Graph generation for low thresholds may take a long time")
        warn.setContentsMargins(0, 0, 0, 0)
        settings_layout.addWidget(warn, 1, 0, 0, -1, Qt.AlignmentFlag.AlignTop)
        settings_layout.addWidget(genes_selection, 2, 0)
        settings_layout.addWidget(genes_path, 2, 1)
        settings_layout.addWidget(expression_selection, 3, 0)
        settings_layout.addWidget(expression_path, 3, 1)
        settings_group.setLayout(settings_layout)

        # Patient selection --------------------------------------------------------------------
        patients_group = QGroupBox("Patient selection")
        patients_layout = QGridLayout()
        patients_group.setLayout(patients_layout)

        self.setWindowTitle("Pathway viewer")

        export_button = QPushButton("Save current view as GEXF")
        export_button.clicked.connect(self.export_gexf)

        right_layout = QVBoxLayout()
        right_layout.addWidget(settings_group)
        right_layout.addWidget(patients_group)
        right_layout.addWidget(export_button)
        right_pane = QWidget()
        right_pane.setMaximumWidth(400)
        right_pane.setLayout(right_layout)

        mainlayout = QGridLayout()
        self.webview = QWebEngineView()
        self.webview.setHtml("""<p style='font-family: sans-serif'>No data loaded.<br>Change the threshold level to trigger graph recreation.</p>""")
        mainlayout.addWidget(self.webview, 0, 0, -1, 4)
        mainlayout.addWidget(right_pane, 0, 4, 1, 1)

        self.setLayout(mainlayout)

    @pyqtSlot()
    def set_genes_db(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import genes database")
        if(not path):
            return
        else:
            self.genes_path = path

        self.ref_coex = make_coexpression_matrix(self.exp_path, self.genes_path)

    @pyqtSlot()
    def set_expression_data(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import genes database")
        if(not path):
            return
        else:
            self.exp_path = path

        self.ref_coex = make_coexpression_matrix(self.exp_path, self.genes_path)

    @pyqtSlot()
    def export_gexf(self):
        if(not self.active_pathway):
            return
        
        path, _ = QFileDialog.getSaveFileName(self, "Export graph information", "", "GEXF graph representation (*.gexf)")
        if(not path):
            return
            
        nx.write_gexf(self.active_pathway.graph, path)

    def svg_changed(self, svg_data):
        self.webview.setContent(svg_data, "image/svg+xml")
    
    def make_graph(self, threshold):
        self.active_pathway = make_pathway_from_thres(threshold, self.ref_coex)
        dot_graph = nx.nx_pydot.to_pydot(self.active_pathway.graph)
        svg = dot_graph.create_svg()
        self.svg_changed(svg)

if __name__ == "__main__":
    app = QApplication([])
    window = Viewer()
    window.show()
    sys.exit(app.exec_())

