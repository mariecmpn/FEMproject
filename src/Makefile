# Makefile Projet fluides incompressibles

# Compiler used
CC = gfortran-11

# Options de compilation
FFLAGS = -O0

# Chemin ou sont stockes les fichiers binaires
ODIR = ../bin
# Chemin ou sera stocke le fichier executable
EXEDIR = ..

# Liste des fichiers binaires
OBJS = $(ODIR)/main.o $(ODIR)/stockage_creux.o $(ODIR)/remplissage_A_L.o $(ODIR)/quadrature_P2.o $(ODIR)/quadrature.o $(ODIR)/recup_node.o $(ODIR)/computation.o $(ODIR)/numerics.o
# Liste des fichiers mod
MODS = stockage_matrice.mod quadrature_P1.mod quadrature_P2.mod to_vtk.mod numerics.mod remplissage_a_l.mod computation.mod

# Librairies used
LIBS = -llapack -lblas

# nom de l'executable
PROG = elements_finis

# regle de compilation pour le fichier executable final
$(PROG) : $(OBJS)
	$(CC) $(FFLAGS) $(OBJS) $(LIBS) -o $(EXEDIR)/$(PROG)

# regle clean pour supprimer les fichiers binaires et executable
clean:
	rm $(ODIR)/*.o

# Regles construction des fichiers binaires
$(ODIR)/main.o: main.f90 stockage_creux.f90 quadrature_P2.f90 quadrature.f90 recup_node.f90 computation.f90 numerics.f90
	$(CC) -c $(FFLAGS) main.f90 -o $(ODIR)/main.o
$(ODIR)/numerics.o: numerics.f90 
	$(CC) -c $(FFLAGS) numerics.f90 -o $(ODIR)/numerics.o
$(ODIR)/stockage_creux.o: stockage_creux.f90 remplissage_A_L.f90 quadrature_P2.f90 quadrature.f90 numerics.f90
	$(CC) -c $(FFLAGS) stockage_creux.f90 -o $(ODIR)/stockage_creux.o
$(ODIR)/quadrature.o: quadrature.f90 numerics.f90
	$(CC) -c $(FFLAGS) quadrature.f90 -o $(ODIR)/quadrature.o
$(ODIR)/recup_node.o: recup_node.f90 numerics.f90
	$(CC) -c $(FFLAGS) recup_node.f90 -o $(ODIR)/recup_node.o
$(ODIR)/quadrature_P2.o: quadrature_P2.f90 quadrature.f90 numerics.f90 
	$(CC) -c $(FFLAGS) quadrature_P2.f90 -o $(ODIR)/quadrature_P2.o
$(ODIR)/remplissage_A_L.o: remplissage_A_L.f90 quadrature_P2.f90 quadrature.f90 numerics.f90
	$(CC) -c $(FFLAGS) remplissage_A_L.f90 -o $(ODIR)/remplissage_A_L.o
$(ODIR)/computation.o: computation.f90 numerics.f90
	$(CC) -c $(FFLAGS) computation.f90 -o $(ODIR)/computation.o