## CCXStressReader

[![latest release](https://img.shields.io/github/release/Mote3D/CCXStressReader.svg)](http://github.com/Mote3D/CCXStressReader/releases/tag/v1.0)

#### Description

CCXStressReader reads element variable output such as stresses and strains from .dat files created by the free finite-element software [CalculiX](http://www.calculix.de/) [G. Dhondt, K. Wittig] or the [FreeCAD](http://www.freecadweb.org/) [J. Riegel, W. Mayer, Y. van Havre] FEM module and computes the minimum, maximum and arithmetic mean of Mises equivalent stress, total effective strain and equivalent plastic strain. The results are written to a .txt file.

The .dat file stores element variable output at the elements' integration points, which generally is more accurate than output extrapolated to the nodes. Evaluating element variable output at the integration points is therefore preferred, particularly in stress analyses with non-linear material behaviour.

Element variable output to the .dat file can be activated by adding the following command to the [CalculiX](http://www.calculix.de/) .inp file:
```
*EL PRINT, ELSET=Eall, FREQUENCY=n
S, E, PEEQ
```
The parameter `FREQUENCY=n` can be used to store element variable output only after every *n*th increment and at the end of the calculation step in non-linear analyses.

#### License

CCXStressReader is licensed under the [GNU Lesser General Public License](https://github.com/Mote3D/CCXStressReader/blob/main/LICENSE.txt).
