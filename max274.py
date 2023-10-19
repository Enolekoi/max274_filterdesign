import sys
import os
import numpy as np
import filterorder as order
import approximation as approx
import widerstand as res 
import plot 
# importiere benötigte Module fuer GUI
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

# Inkludiere Bild in executable
def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

# Caste String zu Float und ersetze leere Strings mit 0.0
def str_to_float(string):
    float_out = float((string[::1].replace(''[::1],'0'[::1], 1))[::1])
    return float_out

# Caste String zu Int und ersetze leere Strings mit 0.0
def str_to_int(string):
    int_out = int((string[::1].replace(''[::1],'0'[::1], 1))[::1])
    return int_out

# Erstelle Hauptfenster Klasse
class Fenster(QWidget):
    # Deklariere Variablen
    Filtertyp   = 'Tiefpass'
    As          = 20        
    Ad          = 3
    fg_tp       = 1000
    fs_tp       = 4000
    fg_hp       = 4000
    fs_hp       = 1000
    dfd         = 1000
    dfs         = 4000
    fM          = 1000
    Nbw         = 2
    Nbe         = 2
    Nch         = 2
    Nca         = 2
    N           = 1
    cnt         = 0
    Approx      = 'Butterworth'
    P           = np.empty(0)
    Z           = np.empty(0)
    f0          = np.empty(0)
    Q           = np.empty(0)
    resP        = np.empty(0)
    resZ        = np.empty(0)
    resf0       = np.empty(0)
    resQ        = np.empty(0)
    f02         = np.empty(0)
    Q2          = np.empty(0)
    FC          = np.empty(0)
    H0          = np.empty(0)
    R1          = np.empty(0)
    R2          = np.empty(0)
    R3          = np.empty(0)
    R4          = np.empty(0)
    
    # Funktionen, die bei Start der Klasse ausgeführt werden
    def __init__(self):
        super().__init__()
        self.initWindow()           # Initialisiere Hauptfenster
        # Ermittle Grundeinstellungen
        self.startup_visibility()   
        self.calc_min_order()
        self.check_export_button() 
        self.check_res_buttons()

    def initWindow(self):

        # Label 
        #Filtertyp
        lb_filter       = QLabel("Filtertyp",self)
        # Daempfung
        lb_daemp    = QLabel("Dämpfung",self)
        lb_Ad       = QLabel("A<sub>d<\sub>",self)
        lb_As       = QLabel("A<sub>s<\sub>",self)
        # Tiefpass
        lb_tp       = QLabel("Tiefpass",self)
        lb_fg_tp    = QLabel("f<sub>g<\sub>",self)
        lb_fs_tp    = QLabel("f<sub>s<\sub>",self)
        # Hochpass
        lb_hp       = QLabel("Hochpass",self)
        lb_fg_hp    = QLabel("f<sub>g<\sub>",self)
        lb_fs_hp    = QLabel("f<sub>s<\sub>",self)
        # Bandpass
        lb_bp       = QLabel("Bandpass",self)
        lb_fM_bp    = QLabel("f<sub>M<\sub>",self)
        lb_dfd_bp   = QLabel("Δf<sub>d<\sub>",self)
        lb_dfs_bp   = QLabel("Δf<sub>s<\sub>",self)

        lb_space    = QLabel("",self)
        # Buttons
        #Filtertyp
        rb_typ_auswahl  = QComboBox(self)
        rb_typ_auswahl.addItem("Tiefpass")
        rb_typ_auswahl.addItem("Hochpass")
        rb_typ_auswahl.addItem("Bandpass")
        rb_typ_auswahl.setToolTip("Auswahl des Filtertypen (nur Tief- und Bandpass auf dem MAX247 umsetzbar)")
        # Filterapproximation 
        rb_f0q_bw       = QRadioButton("Butterworth",self)
        rb_f0q_bw.setToolTip("Butterworth Approximation auswählen")
        rb_f0q_bw.setChecked(True)
        rb_f0q_be       = QRadioButton("Bessel",self)
        rb_f0q_be.setToolTip("Bessel Approximation auswählen")
        rb_f0q_ch       = QRadioButton("Tschebyscheff",self)
        rb_f0q_ch.setToolTip("Tschebyscheff Approximation auswählen")
        rb_f0q_ca       = QRadioButton("Cauer",self)
        rb_f0q_ca.setToolTip("Cauer Approximation auswählen")
        pb_calc_pole    = QPushButton("Filter berechnen",self)
        pb_calc_pole.setToolTip("Berechne Pol- und Nullstellen, sowie" + " f<sub>0</sub> " + "und Q des Filters")
        pb_plot_design  = QPushButton("Plot",self)
        pb_plot_design.setToolTip("Plotte Betrag, Phase und Gruppenlaufzeit der Übertragungsfunktion (Design Sektion)")
        self.pb_export  = QPushButton("Exportieren",self)
        self.pb_export.setToolTip("Exportiere Pol-, Nullstellen, f<sub>0</sub> und Q des Filters in die Sektion zur Berechnung der Widerstandswerte")
        # Line Edit Fields
        # Daempfung
        # self.ef_Ad       = QLineEdit("3",self)
        self.ef_Ad       = QDoubleSpinBox(self)
        self.ef_Ad.setRange(0.0, 100.0)
        self.ef_Ad.setSuffix(" dB")
        self.ef_Ad.setValue(3.0)
        self.ef_Ad.setToolTip("Zugelassener Rippel im Durchlassbereich von Butterworth-, Tschebyscheff- und Cauerfilter (0.01:100.00) [dB] ")
        self.ef_As       = QDoubleSpinBox(self)
        self.ef_As.setRange(0.0, 100.0)
        self.ef_As.setSuffix(" dB")
        self.ef_As.setValue(20.0)
        self.ef_As.setToolTip("Minimal zugelassene Dämpfung im Sperrbereich (0.01:100.00) [dB]")
        # Tiefpass
        self.ef_fg_tp    = QDoubleSpinBox(self)
        self.ef_fg_tp.setRange(1.0, 10000.0)
        self.ef_fg_tp.setSuffix(" Hz")
        self.ef_fg_tp.setValue(1000.0)
        self.ef_fg_tp.setToolTip("Grenzfrequenz des Durchlassbereichs bei Tiefpässen (1:10000) [Hz]")

        self.ef_fs_tp    = QDoubleSpinBox(self)
        self.ef_fs_tp.setRange(1.0, 10000.0)
        self.ef_fs_tp.setSuffix(" Hz")
        self.ef_fs_tp.setValue(4000.0)
        self.ef_fs_tp.setToolTip("Grenzfrequenz des Sperrbereichs bei Tiefpässen (1:10000) [Hz]")
        # Hochpass
        self.ef_fg_hp    = QDoubleSpinBox(self)
        self.ef_fg_hp.setRange(1.0, 10000.0)
        self.ef_fg_hp.setSuffix(" Hz")
        self.ef_fg_hp.setValue(4000.0)
        self.ef_fg_hp.setToolTip("Grenzfrequenz des Durchlassbereichs bei Hochpässen (1:10000) [Hz]")

        self.ef_fs_hp    = QDoubleSpinBox(self)
        self.ef_fs_hp.setRange(1.0, 10000.0)
        self.ef_fs_hp.setSuffix(" Hz")
        self.ef_fs_hp.setValue(1000.0)
        self.ef_fs_hp.setToolTip("Grenzfrequenz des Sperrbereichs bei Hochpässen (1:10000) [Hz]")
        #Bandpass
        self.ef_fM_bp    = QDoubleSpinBox(self)
        self.ef_fM_bp.setRange(1.0, 10000.0)
        self.ef_fM_bp.setSuffix(" Hz")
        self.ef_fM_bp.setValue(1000.0)
        self.ef_fM_bp.setToolTip("Mittenfrequenz eines symetrischen Bandpasses (1:10000) [Hz]")
        self.ef_dfd_bp   = QDoubleSpinBox(self)
        self.ef_dfd_bp.setRange(1.0, 10000.0)
        self.ef_dfd_bp.setSuffix(" Hz")
        self.ef_dfd_bp.setValue(1000.0)
        self.ef_dfd_bp.setToolTip("Bandbreite des Durchlassbereichs eines Bandpasses (1:10000) [Hz]")
        self.ef_dfs_bp   = QDoubleSpinBox(self)
        self.ef_dfs_bp.setRange(1.0, 10000.0)
        self.ef_dfs_bp.setSuffix(" Hz")
        self.ef_dfs_bp.setValue(4000.0)
        self.ef_dfs_bp.setToolTip("Bandbreite des Sperrbereichs eines Bandpasses (1:10000) [Hz]")
        # Filtergrad
        self.ef_n_bw    = QLineEdit(self)
        self.ef_n_bw.setValidator(QIntValidator(0,24,self))
        self.ef_n_bw.setToolTip("Filtergrad eines Butterworth Filters")
        self.ef_n_be    = QLineEdit(self)
        self.ef_n_be.setValidator(QIntValidator(0,24,self))
        self.ef_n_be.setToolTip("Filtergrad eines Bessel Filters")
        self.ef_n_ch    = QLineEdit(self)
        self.ef_n_ch.setValidator(QIntValidator(0,24,self))
        self.ef_n_ch.setToolTip("Filtergrad eines Tschebyscheff Filters")
        self.ef_n_ca    = QLineEdit(self)
        self.ef_n_ca.setValidator(QIntValidator(0,24,self))
        self.ef_n_ca.setToolTip("Filtergrad eines Cauer Filters")
        # Seperator
        sep = QFrame(self)
        sep.setFrameShape(QFrame.HLine)
        sep.setSizePolicy(QSizePolicy.Minimum,QSizePolicy.Expanding)
        sep.setLineWidth(1)
        # Layouts
        g_1     = QGridLayout()
        # Filtertyp
        g_1.addWidget(lb_filter,        1, 2)
        g_1.addWidget(rb_typ_auswahl,   2, 2)
        # Daempfung
        g_1.addWidget(lb_daemp,         3, 2)
        g_1.addWidget(lb_Ad,            4, 1)
        g_1.addWidget(lb_As,            5, 1)
        g_1.addWidget(self.ef_Ad,       4, 2)
        g_1.addWidget(self.ef_As,       5, 2)
        # Tiefpass
        g_1.addWidget(lb_tp,            3, 4)
        g_1.addWidget(lb_fg_tp,         4, 3)
        g_1.addWidget(lb_fs_tp,         5, 3)
        g_1.addWidget(self.ef_fg_tp,    4, 4)
        g_1.addWidget(self.ef_fs_tp,    5, 4)
        # Hochpass
        g_1.addWidget(lb_hp,            3, 6)
        g_1.addWidget(lb_fs_hp,         4, 5)
        g_1.addWidget(lb_fg_hp,         5, 5)
        g_1.addWidget(self.ef_fs_hp,    4, 6)
        g_1.addWidget(self.ef_fg_hp,    5, 6)
        # Bandpass
        g_1.addWidget(lb_bp,            3, 8)
        g_1.addWidget(lb_fM_bp,         4, 7)
        g_1.addWidget(lb_dfd_bp,        5, 7)
        g_1.addWidget(lb_dfs_bp,        6, 7)
        g_1.addWidget(self.ef_fM_bp,    4, 8)
        g_1.addWidget(self.ef_dfd_bp,   5, 8)
        g_1.addWidget(self.ef_dfs_bp,   6, 8)
        # Buttons
        g_1.addWidget(rb_f0q_bw,        8, 2)
        g_1.addWidget(self.ef_n_bw,     9, 2)
        g_1.addWidget(lb_space,         10, 2)
        g_1.addWidget(sep,              11, 2, 1, 8)
        g_1.addWidget(lb_space,         12, 2)
        g_1.addWidget(pb_calc_pole,     13, 2)

        g_1.addWidget(rb_f0q_be,        8, 4)
        g_1.addWidget(self.ef_n_be,     9, 4)
        g_1.addWidget(pb_plot_design,   13, 4)

        g_1.addWidget(rb_f0q_ch,        8, 6)
        g_1.addWidget(self.ef_n_ch,     9, 6)
        g_1.addWidget(self.pb_export,   13, 6)

        g_1.addWidget(rb_f0q_ca,        8, 8) 
        g_1.addWidget(self.ef_n_ca,     9, 8)
    
        # Events, etc.
        # Filtertyp
        rb_typ_auswahl.currentIndexChanged.connect(self.typ_select)
        rb_typ_auswahl.currentIndexChanged.connect(self.input_read)
        rb_typ_auswahl.currentIndexChanged.connect(self.calc_min_order)
        rb_typ_auswahl.currentIndexChanged.connect(self.order_read)
        rb_typ_auswahl.currentIndexChanged.connect(self.calc_f0q)
        rb_typ_auswahl.currentIndexChanged.connect(self.check_export_button)
        # Daempfung
        self.ef_Ad.editingFinished.connect(self.input_read)
        self.ef_As.editingFinished.connect(self.input_read)
        self.ef_As.editingFinished.connect(self.calc_min_order)
        self.ef_Ad.editingFinished.connect(self.calc_min_order)
        #Tiefpass
        self.ef_fg_tp.editingFinished.connect(self.input_read)
        self.ef_fs_tp.editingFinished.connect(self.input_read)
        self.ef_fg_tp.editingFinished.connect(self.calc_min_order)
        self.ef_fs_tp.editingFinished.connect(self.calc_min_order)
        # Hochpass
        self.ef_fg_hp.editingFinished.connect(self.input_read)
        self.ef_fs_hp.editingFinished.connect(self.input_read)
        self.ef_fg_hp.editingFinished.connect(self.calc_min_order)
        self.ef_fs_hp.editingFinished.connect(self.calc_min_order)
        # Bandpass
        self.ef_dfd_bp.editingFinished.connect(self.input_read)
        self.ef_dfs_bp.editingFinished.connect(self.input_read)
        self.ef_fM_bp.editingFinished.connect(self.input_read)
        self.ef_dfs_bp.editingFinished.connect(self.calc_min_order)
        self.ef_dfd_bp.editingFinished.connect(self.calc_min_order)
        self.ef_fM_bp.editingFinished.connect(self.calc_min_order)
        # Filterordnung
        self.ef_n_bw.editingFinished.connect(self.input_read)
        self.ef_n_be.editingFinished.connect(self.input_read)
        self.ef_n_ch.editingFinished.connect(self.input_read)
        self.ef_n_ca.editingFinished.connect(self.input_read)

        pb_calc_pole.pressed.connect(self.order_read)
        pb_calc_pole.pressed.connect(self.calc_f0q)
        pb_calc_pole.pressed.connect(self.check_export_button)
        pb_plot_design.pressed.connect(self.plot_design)
        # Buttons unter Filterordnung
        self.pb_export.pressed.connect(self.export_f0q)
        self.pb_export.pressed.connect(self.calc_res)
        self.pb_export.pressed.connect(self.check_res_buttons)
        rb_f0q_bw.toggled.connect(self.approx_select)
        rb_f0q_be.toggled.connect(self.approx_select)
        rb_f0q_ch.toggled.connect(self.approx_select)
        rb_f0q_ca.toggled.connect(self.approx_select)
        rb_f0q_bw.toggled.connect(self.check_export_button)
        rb_f0q_be.toggled.connect(self.check_export_button)
        rb_f0q_ch.toggled.connect(self.check_export_button)
        rb_f0q_ca.toggled.connect(self.check_export_button)


        # Input Error Feld
        # Label 
        self.lb_inputERROR   = QLabel(self)
        
        # Ordnung Error Feld
        # Label 
        self.lb_orderERROR   = QLabel(self)


        # P-N Feld
        # Label 
        self.lb_pol_null_table  = QLabel('Pol- und Nullstellen', self)
        self.lb_erste_ordnung   = QLabel('Pol erster Ordnung', self)
        self.lb_zweite_ordnung   = QLabel('Pole zweiter Ordnung', self)
        # Table Fields
        self.tb_pol     = QTableWidget(self)
        self.tb_pol.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tb_pol.setColumnCount(2)       
        self.tb_pol.setMaximumHeight(150)
        self.tb_pol.setMaximumWidth(900)
        self.tb_pol.setMinimumHeight(150)
        self.tb_pol.setMinimumWidth(900)
        self.tb_pol.setHorizontalHeaderLabels( ('Polstellen', 'Nullstellen') )
        self.tb_pol.horizontalHeader().setStretchLastSection(True)
        self.tb_pol.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.tb_pol.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)

        self.tb_f0_q_2 = QTableWidget(self)
        self.tb_f0_q_2.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tb_f0_q_2.setColumnCount(2)       
        self.tb_f0_q_2.setMaximumHeight(150)
        self.tb_f0_q_2.setMaximumWidth(450)
        self.tb_f0_q_2.setMinimumHeight(150)
        self.tb_f0_q_2.setMinimumWidth(450)
        self.tb_f0_q_2.setHorizontalHeaderLabels( ("f0", "Q") )
        self.tb_f0_q_2.horizontalHeader().setStretchLastSection(True)
        self.tb_f0_q_2.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.tb_f0_q_2.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)

        self.tb_f0_q_1 = QTableWidget(self)
        self.tb_f0_q_1.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tb_f0_q_1.setColumnCount(2)       
        self.tb_f0_q_1.setMaximumHeight(150)
        self.tb_f0_q_1.setMaximumWidth(450)
        self.tb_f0_q_1.setMinimumHeight(150)
        self.tb_f0_q_1.setMinimumWidth(450)
        self.tb_f0_q_1.setHorizontalHeaderLabels( ("f0", "Q") )
        self.tb_f0_q_1.horizontalHeader().setStretchLastSection(True)
        self.tb_f0_q_1.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.tb_f0_q_1.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)

        v_table_pol         = QVBoxLayout()
        v_table_pol.addWidget(self.lb_pol_null_table)
        v_table_pol.addWidget(self.tb_pol)

        g_table_f0q         = QGridLayout()
        g_table_f0q.addWidget(self.lb_zweite_ordnung, 1, 1)
        g_table_f0q.addWidget(self.lb_erste_ordnung , 1, 2)
        g_table_f0q.addWidget(self.tb_f0_q_2, 2, 1)
        g_table_f0q.addWidget(self.tb_f0_q_1, 2, 2)
        v_table             = QVBoxLayout()
        v_table.addLayout(v_table_pol)
        v_table.addLayout(g_table_f0q)
        # Events, etc.


        ## Widerstandsberechnung
        self.lb_res_Filtertyp = QLabel("",self)
        self.pb_calc_res     = QPushButton("Berechne Widerstandswerte",self)
        self.pb_calc_res.setToolTip("Berechne aus f<sub>0</sub> und Q: \n H<sub>0</sub>, Widerstandswerte in [kΩ]; ermittle mit welchem Potential der Pin FC verbunden werden muss")
        self.pb_calc_pn      = QPushButton("Berechne Pol-/ Nullstellen aus Widerstandswerten",self)
        self.pb_calc_pn.setToolTip("Berechne aus R<sub>1</sub> bis R<sub>4</sub>: f<sub>0</sub>, Q, Pol- und Nullstellen")
        self.pb_plot_res     = QPushButton("Plot",self)
        self.pb_plot_res.setToolTip("Plotte die ausgewählten Sektionen")
        self.cb_plot_sek     = QComboBox(self)
        self.cb_plot_sek.setToolTip("Wähle aus, welche Sektion geplottet werden soll")
        self.cb_plot_sek.addItem('Alle Sektionen')

        # Sektionen + f0 / Q
        self.tb_res_f0_q = QTableWidget(self)
        self.tb_res_f0_q.setColumnCount(8)       
        self.tb_res_f0_q.setMaximumWidth(900)
        self.tb_res_f0_q.setMaximumHeight(150)
        self.tb_res_f0_q.setMinimumWidth(900)
        self.tb_res_f0_q.setMinimumHeight(150)
        self.tb_res_f0_q.setHorizontalHeaderLabels( ("f0", "Q", "H0", "FC verbunden mit", "R1 [kΩ]", "R2 [kΩ]", "R3 [kΩ]", "R4 [kΩ]") )
        self.tb_res_f0_q.horizontalHeader().setStretchLastSection(True)
        self.tb_res_f0_q.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.tb_res_f0_q.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.tb_res_f0_q.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self.tb_res_f0_q.horizontalHeader().setSectionResizeMode(3, QHeaderView.Stretch)
        self.tb_res_f0_q.horizontalHeader().setSectionResizeMode(4, QHeaderView.Stretch)
        self.tb_res_f0_q.horizontalHeader().setSectionResizeMode(5, QHeaderView.Stretch)
        self.tb_res_f0_q.horizontalHeader().setSectionResizeMode(6, QHeaderView.Stretch)
        self.tb_res_f0_q.horizontalHeader().setSectionResizeMode(7, QHeaderView.Stretch)

        self.tb_res_pol     = QTableWidget(self)
        self.tb_res_pol.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tb_res_pol.setColumnCount(2)       
        self.tb_res_pol.setMaximumWidth(900)
        self.tb_res_pol.setMaximumHeight(200)
        self.tb_res_pol.setMinimumWidth(900)
        self.tb_res_pol.setMinimumHeight(200)
        self.tb_res_pol.setHorizontalHeaderLabels( ('Polstellen', 'Nullstellen') )
        self.tb_res_pol.horizontalHeader().setStretchLastSection(True)
        self.tb_res_pol.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.tb_res_pol.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)

        # picture
        schaltung           = QPixmap(resource_path('max_intern.png'))
        schaltung_scaled    = schaltung.scaled(900, 900, Qt.KeepAspectRatio)
        lb_pic_schaltung    = QLabel(self)
        lb_pic_schaltung.setPixmap(schaltung_scaled)
        lb_pic_caption      = QLabel("<i>Abbildung: Beschaltung einer Sektion des MAX247</i>",self)
        # Layout
        g_res_inputs = QGridLayout()
        g_res_inputs.addWidget(self.pb_calc_res, 1, 1)
        g_res_inputs.addWidget(self.pb_calc_pn , 1, 2)
        g_res_inputs.addWidget(self.pb_plot_res, 1, 3)
        g_res_inputs.addWidget(self.cb_plot_sek, 1, 4)
        g_res_inputs.addWidget(self.lb_res_Filtertyp, 2, 1, 1, 4)
        h_res_inputs = QHBoxLayout()
        h_res_inputs.addLayout(g_res_inputs)
        h_res_inputs.addStretch(1)

        h_res_table_f0q = QHBoxLayout()
        h_res_table_f0q.addWidget(self.tb_res_f0_q)
        h_res_table_f0q.addStretch(1)
        h_res_table_pol = QHBoxLayout()
        h_res_table_pol.addWidget(self.tb_res_pol)
        h_res_table_pol.addStretch(1)

        v_res_table = QVBoxLayout()
        v_res_table.addWidget(self.tb_res_f0_q)
        v_res_table.addWidget(self.tb_res_pol)
        v_res_table.addStretch(1)
        # Events, etc.
        self.pb_calc_res.pressed.connect(self.calc_res)
        self.pb_calc_res.pressed.connect(self.check_res_buttons)
        self.pb_calc_pn.pressed.connect(self.calc_f0q_from_r)
        self.pb_calc_pn.pressed.connect(self.check_res_buttons)
        self.pb_plot_res.pressed.connect(self.plot_sek_res)
        # self.pb_plot_res.pressed.connect(self.plot_res)
        self.pb_plot_res.pressed.connect(self.check_res_buttons)

        # Erstelle Design Layout
        tab1_h1= QHBoxLayout()
        tab1_h1.addLayout(g_1)
        tab1_h1.addStretch(1)

        tab1_v2= QVBoxLayout()
        tab1_v2.addWidget(self.lb_inputERROR)

        tab1_v3= QVBoxLayout()
        tab1_v3.addWidget(self.lb_orderERROR)

        tab1_h4= QHBoxLayout()
        tab1_h4.addLayout(v_table)
        tab1_h4.addStretch(1) 
        tab1_main = QVBoxLayout()
        tab1_main.addLayout(tab1_h1)
        tab1_main.addLayout(tab1_v2)
        tab1_main.addLayout(tab1_v3)
        tab1_main.addLayout(tab1_h4)
        tab1_main.addStretch(1)
        
        tab2_h1 = QHBoxLayout()
        tab2_h1.addLayout(v_res_table)
        tab2_h1.addStretch(1)
        tab2_main = QVBoxLayout()
        tab2_main.addLayout(tab2_h1)
        tab2_main.addLayout(h_res_inputs)
        tab2_main.addWidget(lb_pic_schaltung)
        tab2_main.addWidget(lb_pic_caption)
        tab2_main.addStretch(1)
        # Tabs
        self.tabs = QTabWidget(self)
        self.tab1 = QWidget(self)
        self.tab2 = QWidget(self)
        self.tabs.addTab(self.tab1, "Filterdesign")
        self.tabs.addTab(self.tab2, "Widerstandsberechnung")
        
        v_main = QVBoxLayout()
        v_main.addWidget(self.tabs)
        self.tab1.setLayout(tab1_main)
        self.tab2.setLayout(tab2_main)
        
        # Waehle Main Layout fuer Fenster aus
        self.setLayout(v_main)
        # Fenstergeometrie
        self.setGeometry(50, 50, 1200, 1000)
        # Fenstertitel
        self.setWindowTitle("MAX247 Filterdesign")
        # Zeige Fenster
        self.show()
    
    def typ_select(self):
        sender = self.sender()
        self.Filtertyp = sender.currentText()
        if self.Filtertyp == 'Tiefpass':
            #Enable TP Fields
            self.ef_fg_tp.setStyleSheet("background-color:white")
            self.ef_fs_tp.setStyleSheet("background-color:white")
            #Disable HP Fields
            self.ef_fg_hp.setStyleSheet("background-color:lightgrey")
            self.ef_fs_hp.setStyleSheet("background-color:lightgrey")
            #Disable BP Fields
            self.ef_dfs_bp.setStyleSheet("background-color:lightgrey")
            self.ef_dfd_bp.setStyleSheet("background-color:lightgrey")
            self.ef_fM_bp.setStyleSheet("background-color:lightgrey")
        elif self.Filtertyp == 'Hochpass':
            #Disable TP Fields
            self.ef_fg_tp.setStyleSheet("background-color:lightgrey")
            self.ef_fs_tp.setStyleSheet("background-color:lightgrey")
            #Enable HP Fields
            self.ef_fg_hp.setStyleSheet("background-color:white")
            self.ef_fs_hp.setStyleSheet("background-color:white")
            #Disable BP Fields
            self.ef_dfs_bp.setStyleSheet("background-color:lightgrey")
            self.ef_dfd_bp.setStyleSheet("background-color:lightgrey")
            self.ef_fM_bp.setStyleSheet("background-color:lightgrey")
        elif self.Filtertyp == 'Bandpass':
            #Disable TP Fields
            self.ef_fg_tp.setStyleSheet("background-color:lightgrey")
            self.ef_fs_tp.setStyleSheet("background-color:lightgrey")
            #Disable HP Fields
            self.ef_fg_hp.setStyleSheet("background-color:lightgrey")
            self.ef_fs_hp.setStyleSheet("background-color:lightgrey")
            #Enable BP Fields
            self.ef_dfs_bp.setStyleSheet("background-color:white")
            self.ef_dfd_bp.setStyleSheet("background-color:white")
            self.ef_fM_bp.setStyleSheet("background-color:white")

    def approx_select(self):
        sender = self.sender()
        self.Approx = sender.text()
        if self.Approx == 'Butterworth':
            # Enable BW Order Field -- Disable rest
            self.ef_n_bw.setStyleSheet("background-color:white")
            self.ef_n_be.setStyleSheet("background-color:lightgrey")
            self.ef_n_ch.setStyleSheet("background-color:lightgrey")
            self.ef_n_ca.setStyleSheet("background-color:lightgrey")
        elif self.Approx == 'Bessel':
            # Enable BE Order Field -- Disable rest
            self.ef_n_bw.setStyleSheet("background-color:lightgrey")
            self.ef_n_be.setStyleSheet("background-color:white")
            self.ef_n_ch.setStyleSheet("background-color:lightgrey")
            self.ef_n_ca.setStyleSheet("background-color:lightgrey")
        elif self.Approx == 'Tschebyscheff':
            # Enable CH Order Field -- Disable rest
            self.ef_n_bw.setStyleSheet("background-color:lightgrey")
            self.ef_n_be.setStyleSheet("background-color:lightgrey")
            self.ef_n_ch.setStyleSheet("background-color:white")
            self.ef_n_ca.setStyleSheet("background-color:lightgrey")
        elif self.Approx == 'Cauer':
            # Enable CA Order Field -- Disable rest
            self.ef_n_bw.setStyleSheet("background-color:lightgrey")
            self.ef_n_be.setStyleSheet("background-color:lightgrey")
            self.ef_n_ch.setStyleSheet("background-color:lightgrey")
            self.ef_n_ca.setStyleSheet("background-color:white")

    def input_read(self):
            self.As      =self.ef_As.value()
            self.Ad      =self.ef_Ad.value()
            self.fs_tp   =self.ef_fs_tp.value()
            self.fg_tp   =self.ef_fg_tp.value()
            self.fs_hp   =self.ef_fs_hp.value()
            self.fg_hp   =self.ef_fg_hp.value()
            self.dfd     =self.ef_dfd_bp.value()
            self.dfs     =self.ef_dfs_bp.value()
            self.fM      =self.ef_fM_bp.value()

    def order_read(self):
            self.Nbw = str_to_int(self.ef_n_bw.text())
            self.Nbe = str_to_int(self.ef_n_be.text())
            self.Nch = str_to_int(self.ef_n_ch.text())
            self.Nca = str_to_int(self.ef_n_ca.text())

    def calc_min_order(self):
        if (   (self.As    <= self.Ad) 
            or (self.fs_tp <= self.fg_tp) 
            or (self.fg_hp <= self.fs_hp) 
            or (self.dfs   <= self.dfd)
            ):
            self.lb_inputERROR.setText('<b>ERROR: Tolleranzvorgaben fehlerhaft</b>')
            self.ef_n_bw.setText( '0' )
            self.ef_n_be.setText( '0' )
            self.ef_n_ch.setText( '0' ) 
            self.ef_n_ca.setText( '0' )
        else:
            self.lb_inputERROR.setText('')
            if self.Filtertyp == 'Tiefpass':
                N = order.get_lp_filter_order(self.Ad, self.As, self.fg_tp, self.fs_tp)
            elif self.Filtertyp == 'Hochpass':
                N = order.get_hp_filter_order(self.Ad, self.As, self.fg_hp, self.fs_hp)
            elif self.Filtertyp == 'Bandpass':
                N = order.get_bp_filter_order(self.Ad, self.As, self.fg_tp, self.fs_tp)
            # Speicher Ordnung in Variable
            self.Nbw = N[0]
            self.Nbe = N[2]
            self.Nch = N[1]
            self.Nca = N[3]
            # Ueberpruefe Filterordnung
            self.lb_orderERROR.setText('')
            error = False

            if self.Nbw < 2:
                self.Nbw = 2
            elif self.Nbw > 24:
                self.Nbw = 0
                error = True

            if self.Nbe == 1:
                self.Nbe = 2 
            elif self.Nbe > 24:
                self.Nbe = 0
                error = True
            
            if (self.Nch < 2) or (self.Nch > 24):
                self.Nch = 2
            elif self.Nch > 24:
                self.Nch = 0
                error = True

            if self.Nca < 2:
                self.Nca = 2
            elif self.Nca > 24:
                self.Nca = 0
                error = True
            # Bei Fehler -> output Fehlerordnung
            if error == True:
                self.lb_orderERROR.setText('<b>ERROR: ungültige Filterordnung (N &lt; 2 oder N &gt; 24)</b>')
            # setze Filterordnungstext
            self.ef_n_bw.setText( str(self.Nbw) )
            self.ef_n_be.setText( str(self.Nbe) )
            self.ef_n_ch.setText( str(self.Nch) ) 
            self.ef_n_ca.setText( str(self.Nca) )

    def calc_f0q(self):
        if self.Approx == 'Butterworth':
            self.N = int(self.Nbw)
        elif self.Approx == 'Bessel':
            self.N = int(self.Nbe)
        elif self.Approx == 'Tschebyscheff':
            self.N = int(self.Nch)
        elif self.Approx == 'Cauer':
            self.N = int(self.Nca)
        if (self.N < 2) or (self.N > 24):
            self.lb_orderERROR.setText('<b>Error: ungültige Filterordnung (N &lt; 2 oder N &gt; 24); Berechnung nicht möglich</bf>')
        else:
            self.lb_orderERROR.setText('')
            if self.Filtertyp == 'Tiefpass':
                self.P,self.Z,self.f0,self.Q =  \
                    approx.lp_approximations(self.Ad, self.As, self.fg_tp, self.fs_tp,self.N,self.Approx)
            elif self.Filtertyp == 'Hochpass':
                self.P,self.Z,self.f0,self.Q =  \
                    approx.hp_approximations(self.Ad, self.As, self.fg_hp, self.fs_hp,self.N,self.Approx)
            elif self.Filtertyp == 'Bandpass':
                self.P,self.Z,self.f0,self.Q =  \
                    approx.bp_approximations(self.Ad, self.As, self.dfd, self.dfs,self.fM,self.N,self.Approx)
            f01, Q1, self.f02, self.Q2 = approx.get_second_first_order(self.f0, self.Q, self.N)
            # anzahl an Elementen in Tabelle
            self.tb_pol.setRowCount(self.N)
            self.tb_f0_q_2.setRowCount(self.f02.size)
            # Ausgabe 
            for n in range(self.P.size):
                self.tb_pol.setItem(n, 0, QTableWidgetItem(str("{:.4f}".format(self.P[n]) ) ) )
            for n in range(self.Z.size):
                self.tb_pol.setItem(n, 1, QTableWidgetItem(str("{:.4f}".format(self.Z[n]) ) ) )
            for n in range(self.f02.size):
                self.tb_f0_q_2.setItem(n, 0, QTableWidgetItem(str("{:.4f}".format(self.f02[n]) ) ) )
                self.tb_f0_q_2.setItem(n, 1, QTableWidgetItem(str("{:.4f}".format(self.Q2[n]) ) ) )
            if len(f01) != 0:
                self.tb_f0_q_1.setRowCount(1)       
                self.tb_f0_q_1.setItem(0, 0, QTableWidgetItem(str("{:.4f}".format(f01[0]) ) ) )
                self.tb_f0_q_1.setItem(0, 1, QTableWidgetItem(str("{:.4f}".format(Q1[0]) ) ) )
            else: 
                self.tb_f0_q_1.setRowCount(0)       

    def check_export_button(self):
        # Ueberpruefe, ob export button gezeigt wird
        if (  (self.tb_f0_q_1.rowCount() > 0)
            or(self.Filtertyp == 'Hochpass' )
            or(self.Approx == 'Cauer')
            or(self.N < 2)
            or(self.N > 24)
            ):
            self.pb_export.setHidden(True)
        else:
            self.pb_export.setHidden(False)

    def export_f0q(self):
        #Exportiere f0, Q, P, Z in Widerstandswert Sektion
        self.lb_res_Filtertyp.setText("Filtertyp: "+self.Filtertyp)
        self.tb_res_f0_q.setRowCount(self.f02.size)       
        self.tb_res_pol.setRowCount(self.P.size) 
        self.resf0 = self.f02
        self.resQ  = self.Q2
        self.resP = self.P
        self.resZ = self.Z
        # Schreibe in Tabellen
        for n in range(self.f02.size):
            self.tb_res_f0_q.setItem(n, 0, QTableWidgetItem(str("{:.4f}".format(self.f02[n]) ) ) )
            self.tb_res_f0_q.setItem(n, 1, QTableWidgetItem(str("{:.4f}".format(self.Q2[n]) ) ) )
            self.tb_res_f0_q.setItem(n, 2, QTableWidgetItem(  ) )
            self.tb_res_f0_q.setItem(n, 3, QTableWidgetItem(  ) )
            self.tb_res_f0_q.setItem(n, 4, QTableWidgetItem(  ) )
            self.tb_res_f0_q.setItem(n, 5, QTableWidgetItem(  ) )
            self.tb_res_f0_q.setItem(n, 6, QTableWidgetItem(  ) )
            self.tb_res_f0_q.setItem(n, 7, QTableWidgetItem(  ) )
        for n in range(self.P.size):
            self.tb_res_pol.setItem(n, 0, QTableWidgetItem(str("{:.4f}".format(self.P[n]) ) ) )
        for n in range(self.Z.size):
            self.tb_res_pol.setItem(n, 1, QTableWidgetItem(str("{:.4f}".format(self.Z[n]) ) ) )
    def calc_res(self):
        # Berechne Widerstandswerte aus f0, Q
        self.R1, self.R2, self.R3, self.R4, self.H0, self.FC =res.berechne_widerstand(self.f02.size
                                                                                      , self.Filtertyp
                                                                                      , self.f02
                                                                                      , self.Q2
                                                                                      , 1)
        # Schreibe in Tabelle
        for n in range(self.f02.size):
            self.tb_res_f0_q.setItem(n, 2, QTableWidgetItem(str("{:.4f}".format(self.H0[n]) ) ) )
            self.tb_res_f0_q.setItem(n, 3, QTableWidgetItem(str(self.FC[n]) ) )
            self.tb_res_f0_q.setItem(n, 4, QTableWidgetItem(str("{:.4f}".format(self.R1[n]) ) ) )
            self.tb_res_f0_q.setItem(n, 5, QTableWidgetItem(str("{:.4f}".format(self.R2[n]) ) ) )
            self.tb_res_f0_q.setItem(n, 6, QTableWidgetItem(str("{:.4f}".format(self.R3[n]) ) ) )
            self.tb_res_f0_q.setItem(n, 7, QTableWidgetItem(str("{:.4f}".format(self.R4[n]) ) ) )

    def calc_f0q_from_r(self):
        # Lese Wiederstandswerte aus Tabelle ein
        for n in range(self.R1.size):
            self.R1[n] = str_to_float(self.tb_res_f0_q.item(n,4).text() )
            self.R2[n] = str_to_float(self.tb_res_f0_q.item(n,5).text() )
            self.R3[n] = str_to_float(self.tb_res_f0_q.item(n,6).text() )
            self.R4[n] = str_to_float(self.tb_res_f0_q.item(n,7).text() )
        # Berechne f0, Q, P, Z aus Widerstandswerten
        self.resf02, self.resQ2, self.resP, self.resZ = res.berechne_pole( self.R1
                                                                          ,self.R2
                                                                          ,self.R3
                                                                          ,self.R4
                                                                          ,self.FC
                                                                          ,self.Filtertyp)
        # Schreibe in Tabelle
        for n in range(self.resf02.size):
            self.tb_res_f0_q.setItem(n, 0, QTableWidgetItem(str("{:.4f}".format(self.resf02[n]) ) ) )
            self.tb_res_f0_q.setItem(n, 1, QTableWidgetItem(str("{:.4f}".format(self.resQ2[n]) ) ) )
        for n in range(self.resP.size):
            self.tb_res_pol.setItem(n, 0, QTableWidgetItem(str("{:.4f}".format(self.resP[n]) ) ) )
        for n in range(self.resZ.size):
            self.tb_res_pol.setItem(n, 1, QTableWidgetItem(str("{:.4f}".format(self.resZ[n]) ) ) )

    def plot_design(self):
        title = 'Übertragungsfunktion ' + self.Approx + ' ' + self.Filtertyp + ' (Design)'
        plot.plot(self.P, self.Z, title)

    def plot_sek_res(self):
        # Alle Sektionen sollen geplottet werden
        if self.cb_plot_sek.currentText() == 'Alle Sektionen':
            title = 'Übertragungsfunktion ' + self.Approx + ' ' + str(self.Filtertyp) + ' (Widerstandsberechnung)'
            plot.plot(self.resP, self.resZ, title)
        else:
            #Bestimmte Sektion soll geplottet werden
            Sek = str_to_int(self.cb_plot_sek.currentText() )
            Sek = Sek-1
            sekP, sekZ = res.f0Q_to_PZ(self.resf0[Sek], self.resQ[Sek], self.Filtertyp)
            title = 'Übertragungsfunktion ' + self.Approx + ' ' + self.Filtertyp + ' (Sektion ' + str(Sek+1) + ')'
            plot.plot(sekP, sekZ, title)
    
    # Plotte die Übertragungsfunktion aus Tab "Widerstandsberechnungs"
    def plot_res(self):
        title = 'Übertragungsfunktion ' + self.Approx + ' ' + self.Filtertyp + ' (Widerstandsberechnung)'
        plot.plot(self.resP, self.resZ, title)

    # Überprüfe welche der Buttons im Tab "Widerstandsberechnung" sichtbar sind
    def check_res_buttons(self):
        self.pb_calc_pn.setHidden(True)
        self.pb_plot_res.setHidden(True)
        self.cb_plot_sek.setHidden(True)
        self.pb_calc_res.setHidden(True)
        
        if self.R1.size > 0:
            # mache "Berechne Pol-/Nullstellen aus Widerstandswerten" sichtbar
            self.pb_calc_pn.setHidden(False)

        if self.tb_res_f0_q.rowCount() > 0:
            # mache "Berechne Widerstandswerte" sichtbar
            self.pb_calc_res.setHidden(False)

        if self.tb_res_pol.rowCount() > 0:
            # mache "Plot" sichtbar
            self.pb_plot_res.setHidden(False)
            self.cb_plot_sek.setHidden(False)

        self.cb_plot_sek.clear()
        self.cb_plot_sek.addItem('Alle Sektionen')

        for n in range(self.tb_res_f0_q.rowCount()):
            self.cb_plot_sek.addItem(str(n+1))

    def startup_visibility(self):
        # Enable BW Order Field -- Disable rest
        self.ef_n_bw.setStyleSheet("background-color:white")
        self.ef_n_be.setStyleSheet("background-color:lightgrey")
        self.ef_n_ch.setStyleSheet("background-color:lightgrey")
        self.ef_n_ca.setStyleSheet("background-color:lightgrey")
        # Enable TP Fields
        self.ef_fg_tp.setStyleSheet("background-color:white")
        self.ef_fs_tp.setStyleSheet("background-color:white")
        # Disable HP Fields
        self.ef_fg_hp.setStyleSheet("background-color:lightgrey")
        self.ef_fs_hp.setStyleSheet("background-color:lightgrey")
        # Disable BP Fields
        self.ef_dfs_bp.setStyleSheet("background-color:lightgrey")
        self.ef_dfd_bp.setStyleSheet("background-color:lightgrey")
        self.ef_fM_bp.setStyleSheet("background-color:lightgrey")

# Starte Application mit Übergabe Parameter
app = QApplication(sys.argv)

w = Fenster()

# Beende Programm, wenn Application geschlossen wird.
sys.exit( app.exec_() )
