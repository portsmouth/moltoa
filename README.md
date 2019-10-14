# moltoa
Python script to convert molecular structures in SDF format to Arnold renderer

Usage is:
```
python moltoa.py my_molecule.sdf
python moltoa.py --style toon my_molecule.sdf
```
This generates `my_molecule.ass` in the working directory, which can be rendered in Arnold with kick via e.g.:
```
kick -nostdin -device cpu -as 4 -ipr -i my_molecule.ass
```

![](./renders/THC.png?raw=true "THC (toon render)")

<img src="./renders/THC.png"/><img src="./renders/THC_toon.png"/>

<img src="./renders/ibuprofen.png"/><img src="./renders/ibuprofen_toon.png"/>

<img src="./renders/caffeine.png"/><img src="./renders/caffeine_toon.png"/>

<img src="./renders/fluoxetine.png"/><img src="./renders/fluoxetine_toon.png"/>
