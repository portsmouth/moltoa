# moltoa
Python script to convert molecular structures in SDF format (available e.g. at https://pubchem.ncbi.nlm.nih.gov) to Arnold renderer

Usage is:
```
python moltoa.py my_molecule.sdf
python moltoa.py --style toon my_molecule.sdf
```
This generates `my_molecule.ass` in the working directory, which can be rendered in Arnold with kick via e.g.:
```
kick -nostdin -device cpu -as 4 -ipr -i my_molecule.ass
```

<img src="./renders/THC.png?raw=true" alt="THC" width="400" height="400"><img src="./renders/THC_toon.png?raw=true" alt="THC toon" width="400" height="400">

<img src="./renders/ibuprofen.png?raw=true" alt="ibuprofen" width="400" height="400"><img src="./renders/ibuprofen_toon.png?raw=true" alt="ibuprofen toon" width="400" height="400">

<img src="./renders/caffeine.png?raw=true" alt="caffeine" width="400" height="400"><img src="./renders/caffeine_toon.png?raw=true" alt="caffeine toon" width="400" height="400">

<img src="./renders/fluoxetine.png?raw=true" alt="fluoxetine" width="400" height="400"><img src="./renders/fluoxetine_toon.png?raw=true" alt="fluoxetine toon" width="400" height="400">

