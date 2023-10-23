# MAX274 Filterdesign Software

## Installation (Linux)
  1. Installiere Abhängigkeiten:

  ### Unter Arch Linux:
  ```
  $ pacman -S python3 qt5 pyqt5 numpy scipy matplotlib pyinstall git
  ```
  ### Unter Debian/Ubuntu:
  ```
  $ apt-get install python3 python3-pip qt5 git
  $ python -m pip install pyqt5 numpy scipy matplotlib pyinstall
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
  $ python -m PyInstaller --onefile --windowed --add-data "max_intern.png:." max274.py
  ```
  5. Die Binärdatei liegt im Verzeichnis 'max274_filterdesign/dist"
  
## Installation (Windows)
  1. Installiere Abhängigkeiten:
  [python3](https://www.python.org/downloads/windows/), [pip](https://pip.pypa.io/en/stable/installation/), [qt5](https://doc.qt.io/qt-5/windows.html), [git](https://git-scm.com/download/win)

  ```
  py -m pip install pyqt5, numpy, scipy, matplotlib, pyinstall
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
  $ python -m PyInstaller --onefile --windowed --add-data "max_intern.png;." max274.py
  ```
  5. Die Binärdatei liegt im Verzeichnis 'max274_filterdesign/dist"
