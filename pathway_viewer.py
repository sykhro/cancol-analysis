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


class PathwayView(QWebEngineView):
    def __init__(self, pathway):
        super().__init__()
        if pathway:
            self.active_pathway = pathway
            self.setContent(nx.nx_pydot.to_pydot(self.active_pathway.graph).create_svg(), "image/svg+xml")
        else:
            self.active_pathway = None
            self.setHtml(
                """<p style='font-family: sans-serif'>No data loaded.<br>Change the threshold level to trigger graph recreation.</p>"""
            )

    def svg_changed(self, svg_data):
        self.setContent(svg_data, "image/svg+xml")


class Viewer(QWidget):
    def __init__(self):
        super().__init__()

        self.exp_path = DEFAULT_EXP
        self.genes_path = DEFAULT_GENES
        self.ref_coex = make_coexpression_matrix(self.exp_path, self.genes_path)

        # Settings -----------------------------------------------------------------------------
        settings_group = QGroupBox("Pathway generation settings")
        thrs_layout = QHBoxLayout()
        self.thrs_box = QDoubleSpinBox()
        self.thrs_box.setMinimum(0)
        self.thrs_box.setMaximum(1.0)
        self.thrs_box.setSingleStep(0.05)
        self.thrs_box.valueChanged.connect(self.make_graph)
        thrs_label = QLabel("Link threshold")
        thrs_layout.addWidget(thrs_label)
        thrs_layout.addWidget(self.thrs_box)

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
        settings_layout.addWidget(self.thrs_box, 0, 1)
        warn = QLabel(
            "Warning! Graph generation for low thresholds may take a long time"
        )
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
        patients_layout.addWidget(QLabel("Patient ID"), 0, 0)
        
        self.patients_dropdown = QComboBox()
        self.patients_dropdown.setEditable(True)
        self.init_patients_list()
        self.add_patient_button = QPushButton("Add patient")
        self.add_patient_button.clicked.connect(self.add_patient)
        self.observed_patients = QListWidget()
        
        patients_layout.addWidget(self.patients_dropdown, 0, 1)
        patients_layout.addWidget(self.add_patient_button, 1, 0, 1, -1)
        patients_layout.addWidget(self.observed_patients, 2, 0, 1, -1)
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
        self.tabs = QTabWidget()
        self.tabs.addTab(
            PathwayView(None),
            "Baseline {}".format(self.thrs_box.value() if self.thrs_box.value() != 0 else ""),
        )
        mainlayout.addWidget(self.tabs, 0, 0, -1, 4)
        mainlayout.addWidget(right_pane, 0, 4, 1, 1)

        self.setLayout(mainlayout)

    def init_patients_list(self):
        self.patients_log = pd.read_csv("TRIBE2_db.csv")
        self.patients_log.sort_values("dpfs", inplace=True)
        patients_list = self.patients_log["PatientFirstName"].to_list()
        model = QStringListModel(patients_list)
        self.patients_dropdown.setModel(model)

    def add_patient(self):
        new_id, new_idx = [self.patients_dropdown.currentText(), self.patients_dropdown.currentIndex()]
        if not self.observed_patients.findItems(new_id, Qt.MatchFlag.MatchExactly):
            pathway = make_pathway_from_thres(self.thrs_box.value(), self.ref_coex)
            new_view = PathwayView(pathway)
            self.tabs.addTab(new_view, pathway.name)
            self.observed_patients.addItem(new_id)
            self.patients_dropdown.removeItem(new_idx)

    @pyqtSlot()
    def set_genes_db(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import genes database")
        if not path:
            return
        else:
            self.genes_path = path

        self.ref_coex = make_coexpression_matrix(self.exp_path, self.genes_path)

    @pyqtSlot()
    def set_expression_data(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import genes database")
        if not path:
            return
        else:
            self.exp_path = path

        self.ref_coex = make_coexpression_matrix(self.exp_path, self.genes_path)

    @pyqtSlot()
    def export_gexf(self):
        pw_view = self.tabs.currentWidget()
        active_pathway = pw_view.active_pathway
        if not active_pathway:
            return

        path, _ = QFileDialog.getSaveFileName(
            self, "Export graph information", "", "GEXF graph representation (*.gexf)"
        )
        if not path:
            return

        nx.write_gexf(active_pathway.graph, path)

    def make_graph(self, threshold):
        pw_view = self.tabs.currentWidget()
        pw_view.active_pathway = make_pathway_from_thres(threshold, self.ref_coex)
        dot_graph = nx.nx_pydot.to_pydot(pw_view.active_pathway.graph)
        pw_view.svg_changed(dot_graph.create_svg())


if __name__ == "__main__":
    app = QApplication([])
    window = Viewer()
    window.show()
    sys.exit(app.exec_())
