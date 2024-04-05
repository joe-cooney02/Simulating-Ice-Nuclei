Files and what they do:
main.py is the main file for the simulation. This includes the summoning of molecules, I/O, and the MMC algorithm.
interactions.py contains the functions used to calculate the energies involved in the system.
cartPolar.py is a helper module that contains functions for coordinate transformations.
cone.py and oval.py define the water molecules and the nucleus, respectively.

mainMT.py is a version of the main.py file that uses multithreading.
main-Joes-Desktop.py is a backup of the main.py file.

Input:
Change the values at the top of main.py to adjust simulation parameters:
    TEMP: controls temperature of system.
        (be aware of the difference in phase diagram between TIP4P-2005 and experimental).
    numMols: the number of water molecules to be simulated.
    spacing: useful for different geometries such as cubic lattices.
    max_steps: number of simulation steps
    nSize: size of nucleus (in angstroms)
    nCharge: charge of the nucleus (in units of elementary charge)

Feel free to adjust other parameters as well (initial geometry, timestep, crystallization condition, etc.)

Output:
When you run main.py, a browser window will open. This is due to the vpython dependency in cone.py and oval.py.
The window will contain a visual representation of the system, text will be output to the console confirming connection.
Four .txt files will be created (Nearest-Neighbor angle, Nearest-Neighbor dist, Number crystallized, Nucleus distance).
These files will contain the extracted values of the system corresponding to their names.

Dependencies, External modules, and Conventions:
You need to install numpy and vpython, math and time are built-in modules.
VPython requires connection to a browser for the simulation to run, this can be accomplished via port-forwarding.
The coordinate system follows the physics standard (phi is azimuthal angle, theta polar).
The unit system follows the Quantum Mechanics conventions (Electron Volts, Angstroms, Elementary Charges)
