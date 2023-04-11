# MADE Plugin
<img src=".MADE_plugin/UI/Icons/MADE_logo.png" width="200">



## Macromolecular (Maribor) Density and Structure Analysis plugin



## Description
MADE plugin is a PyMOL GUI tool that applies the multiple different superposition algorithm for identification and prediction of non-protein species in macromolecular complexes. We have show that the plugin is capable of identifiyng metal ions in their binding sites in *apo* protein structures.

It is based upon the ProBiS H2O plugin (http://insilab.org/probis-h2o/) for detecting conserved water molecules.

## Installation

The MADE plugin is suppported on the Windows and Linux operating systems.

MADE plugin requires the PyMOL (https://pymol.org/) molecular visualization system, the Open-Source PyMOL is available free of charge. Installation instructions can be found on the PyMOL wiki: 

https://pymolwiki.org/index.php/Windows_Install

https://pymolwiki.org/index.php/Linux_Install

MADE plugin also requires the scikit-learn python library (https://scikit-learn.org/), it can be installed though pythons' pip:

    pip install scikit-learn

With PyMOL and scikit-learn installed, install the MADE plugin in the plugins tab of PyMOL:
- Navigate to **Plugin -> Plugin Manager -> Install New Plugin -> Install from local file** and select the **MADE_Plugin.py** file
- Upon first startup of the plugin input location of the **.MADE_Plugin** directory, 
- Selecet a **Local Database Directory** and press **Setup database** in the Settings tab

For a more detailed explanation of the installation and usage of the MADE plugin consult the tutorial pdf.


To reinstall the plugin either delete the .MADE_plugin directory or the .MADE_plugin_installdir.txt file in your home folder and launch the plugin again.

## Citation
Please cite the following paper if you use the MADE plugin in your work:


