# MAX274 Filterdesign Software

## Installation (Linux)
  1. Installiere Abhängigkeiten:
  ```
  python3, qt5, pyqt5, numpy, scipy, matplotlib, pyinstall
  ```
  2. Klone Repository
  ```
  $ git clone https://github.com/Enolekoi/max274_filterdesign.git
  ```
  3. Wechsle ins Arbeitsverzeichnis
  ```
  $ cd /max274_filterdesign
  ```
  4. Erstelle executable mit folgendem Kommando
  ```
  $ python -m PyInstaller --onefile --windowedd --add-data "max_intern.png:." max274.py
  ```
  5. Die Binärdatei liegt im Verzeichnis 'max274_filterdesign/dist"
  
## Installation (Windows)
  1. Installiere Abhängigkeiten:
  ```
  python3, qt5, pyqt5, numpy, scipy, matplotlib, pyinstall
  ```
  2. Wechsle ins Arbeitsverzeichnis
  ```
  $ cd /max274_filterdesign
  ``` 
  3. Erstelle executable mit folgendem Kommando
  ```
  $ python -m PyInstaller --onefile --windowedd --add-data "max_intern.png;." max274.py
  ```
  4. Die Binärdatei liegt im Verzeichnis 'max274_filterdesign/dist"
