#!/usr/bin/env python

import sys

import networkx as nx
import pandas as pd
from PyQt5.QtCore import QStringListModel, Qt, pyqtSlot
from PyQt5.QtGui import QColor
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import (
    QApplication,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from analysis_nx import retrieve_mutations
from network_builder import Pathway, make_pathway_from_thres

DEFAULT_EXP = "GSE40367-stripped.txt"
DEFAULT_GENES = "GPL570-stripped.txt"


def lerp(a, b, t):
    return a * (1 - t) + b * t


def percentage_to_rgb(percentage, base_color=0xFFFFFF, target_color=0xFF0000):
    start = QColor(base_color).toHsv()
    end = QColor(target_color).toHsv()

    h = lerp(start.hue(), end.hue(), percentage / 100)
    s = lerp(start.saturation(), end.saturation(), percentage / 100)
    v = lerp(start.value(), end.value(), percentage / 100)
    rgb = QColor.fromHsv(int(h), int(s), int(v)).getRgb()

    return "#%02x%02x%02x" % (rgb[0], rgb[1], rgb[2])


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
    def __init__(self, pathway, mutations):
        super().__init__()
        self.mutations = mutations
        if pathway:
            self.active_pathway = pathway
            self.setContent(
                nx.nx_pydot.to_pydot(self.active_pathway.graph).create_svg(),
                "image/svg+xml",
            )
        else:
            self.active_pathway = None
            self.setHtml(
                """<p style='font-family: sans-serif'>No data loaded.<br>Set a threshold and click 'Refresh' to trigger graph recreation.</p>"""
            )

    def refresh_svg(self):
        self.setContent(
            nx.nx_pydot.to_pydot(self.active_pathway.graph).create_svg(),
            "image/svg+xml",
        )

    def svg_changed(self, svg_data):
        self.setContent(svg_data, "image/svg+xml")


class Viewer(QWidget):
    def __init__(self):
        super().__init__()

        self.exp_path = DEFAULT_EXP
        self.genes_path = DEFAULT_GENES
        self.sequencing_data = pd.read_csv("TRIBE2_seq_res.csv")
        self.ref_coex = make_coexpression_matrix(self.exp_path, self.genes_path)

        # Settings -----------------------------------------------------------------------------
        settings_group = QGroupBox("Pathway generation settings")
        thrs_layout = QHBoxLayout()
        self.thrs_box = QDoubleSpinBox()
        self.thrs_box.setMinimum(0)
        self.thrs_box.setMaximum(1.0)
        self.thrs_box.setSingleStep(0.05)
        self.thrs_box.setValue(0.6)
        refresh_button = QPushButton("Refresh")
        refresh_button.clicked.connect(self.refresh_view)
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
        settings_layout.addWidget(refresh_button, 0, 2)
        warn = QLabel(
            "Warning! Graph generation for low thresholds may take a long time"
        )
        settings_layout.addWidget(warn, 1, 0, 1, -1, Qt.AlignmentFlag.AlignTop)
        settings_layout.addWidget(genes_selection, 2, 0)
        settings_layout.addWidget(genes_path, 2, 1, 1, 2)
        settings_layout.addWidget(expression_selection, 3, 0)
        settings_layout.addWidget(expression_path, 3, 1, 1, 2)
        settings_group.setLayout(settings_layout)

        # Patient selection --------------------------------------------------------------------
        patients_group = QGroupBox("Patient selection")
        patients_layout = QGridLayout()
        patients_layout.addWidget(QLabel("Patient ID"), 0, 0)

        self.patients_dropdown = QComboBox()
        self.patients_dropdown.setEditable(True)
        self.init_patients_list()
        self.add_patient_button = QPushButton("Add patient ↓")
        self.add_patient_button.clicked.connect(self.add_patient)
        self.remove_patient_button = QPushButton("Remove patient ↑")
        self.remove_patient_button.clicked.connect(self.remove_patient)
        self.observed_patients = QListWidget()

        patients_layout.addWidget(self.patients_dropdown, 0, 1)
        patients_layout.addWidget(self.add_patient_button, 1, 0)
        patients_layout.addWidget(self.remove_patient_button, 1, 1)
        patients_layout.addWidget(self.observed_patients, 3, 0, 1, -1)
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
        self.tabs.tabCloseRequested.connect(lambda index: self.tabs.removeTab(index))
        baseline_idx = self.tabs.addTab(
            PathwayView(None, []),
            "Baseline".format(
                self.thrs_box.value() if self.thrs_box.value() != 0 else ""
            ),
        )
        self.tab_indices = {"Baseline": baseline_idx}
        mainlayout.addWidget(self.tabs, 0, 0, -1, 4)
        mainlayout.addWidget(right_pane, 0, 4, 1, 1)

        self.setLayout(mainlayout)

    def init_patients_list(self):
        self.patients_log = pd.read_csv("TRIBE2_db.csv")
        self.patients_log.sort_values("dpfs", inplace=True)
        patients_list = self.patients_log["PatientFirstName"].to_list()
        model = QStringListModel(patients_list)
        self.patients_dropdown.setModel(model)

    def remove_patient(self):
        idx_remove = self.observed_patients.currentIndex()
        to_close = self.observed_patients.itemAt(
            idx_remove.column(), idx_remove.row()
        ).data(Qt.DisplayRole)
        self.observed_patients.takeItem(idx_remove.row())

        self.tabs.removeTab(self.tab_indices[to_close])
        del self.tab_indices[to_close]

    def add_patient(self):
        new_id = self.patients_dropdown.currentText()
        new_idx = self.patients_dropdown.currentIndex()

        if not self.observed_patients.findItems(new_id, Qt.MatchFlag.MatchExactly):
            mutations = retrieve_mutations(new_id, self.sequencing_data)
            pathway = self.make_annotated_pathway(mutations)

            new_view = PathwayView(pathway, mutations)
            tab_idx = self.tabs.addTab(new_view, new_id)
            self.tab_indices[new_id] = tab_idx
            self.tabs.setCurrentIndex(tab_idx)
            self.observed_patients.addItem(new_id)
            self.patients_dropdown.setCurrentIndex(new_idx + 1)

    def make_annotated_pathway(self, mutations):
        thres = self.thrs_box.value()
        pathway = make_pathway_from_thres(thres, self.ref_coex)
        muts = mutations.copy()

        colors = {}
        if "NGS_PercentMutated" in muts:
            muts["NGS_PercentMutated"] = muts["NGS_PercentMutated"].map(
                percentage_to_rgb
            )
            colors = muts.set_index("Biomarker").to_dict()["NGS_PercentMutated"]

        for node in pathway.graph.nodes:
            pathway.graph.nodes[node]["color"] = "black"
            pathway.graph.nodes[node]["style"] = "filled"

            if node in colors:
                pathway.graph.nodes[node]["fillcolor"] = colors[node]
            else:
                pathway.graph.nodes[node]["fillcolor"] = "#ffffff"

        return pathway

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
            self,
            "Export graph information",
            self.tabs.tabText(self.tabs.currentIndex()),
            "GEXF graph representation (*.gexf)",
        )
        if not path:
            return

        nx.write_gexf(active_pathway.graph, path)

    def refresh_view(self):
        pw_view = self.tabs.currentWidget()
        pw_view.active_pathway = self.make_annotated_pathway(pw_view.mutations)
        pw_view.refresh_svg()


if __name__ == "__main__":
    app = QApplication([])
    window = Viewer()
    window.show()
    sys.exit(app.exec_())
