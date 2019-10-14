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
![](./renders/THC.png?raw=true "THC") ![](./renders/THC_toon.png?raw=true "THC (toon render)")

![](./renders/ibuprofen.png?raw=true "ibuprofen") ![](./renders/ibuprofen_toon.png?raw=true "ibuprofen (toon render)")

![](./renders/caffeine.png?raw=true "caffeine") ![](./renders/caffeine_toon.png?raw=true "caffeine (toon render)")

![](./renders/fluoxetine.png?raw=true "fluoxetine") ![](./renders/fluoxetine_toon.png?raw=true "fluoxetine (toon render)")

