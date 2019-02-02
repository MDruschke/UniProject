# UniProject
In diesem Repository sind 2 repräsentative Beispiele meiner Programmierarbeit, welche wie
in meiner Bewerbung beschrieben, sich größtenteils auf die Datenanalyse großer
Simulation-Outputs und auf die Bildausgabe spezialisiert.
Der Code ist in Python 2.7 und im Style guide PEP 8 geschrieben. Zum besseren Verständnis
habe ich noch ein paar zusätzliche Kommentare eingefügt.
Die Datei "Plot_interpolation.py" ist ein Beispiel der Bilverarbeitung.
In dieser Datei erstelle ich eine Plot, welcher vergleichbar mit Figure 3 & Figure 4 aus
meinem kürzlich veröffentlichten Paper (https://arxiv.org/pdf/1809.10038.pdf) ist. Darin 
wird ein sogenanter "Halo" (ein riesige Ansammlung an Gas und Dunkler Materie im frühen 
Universum, aus welcher sich Sterne bilden können) geplottet. Dabei wird die Dichte farblich 
und die Geschwindigkeiten mithilfe von Pfeilen dargestellt.
Die Datei "Spin_calculation.py" ist ein Beispiel für die Datenanalyse großer 
Simulation-Outputs. In dieser berechne für über 9000 Haloes aus den Gas- und dunkle Materie
Teilchen deren Rotation (Spin), welche ein wichtiges Maß für die Entstehung eines Sternes ist.
Dazu werden u.a. die Position, Masse und Geschwindigkeit jedes Teilchens verwendet. 
In meinem Paper ist die Verteilung dieserSpins in Figure 1 dargestellt.

