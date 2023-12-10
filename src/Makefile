# Makefile Projet fluides incompressibles

# Compiler used
CC = gfortran

# Options de compilation
FFLAGS = -O0

# Chemin ou sont stockes les fichiers binaires
ODIR = ../bin
# Chemin ou sera stocke le fichier executable
EXEDIR = ..

# Liste des fichiers binaires
OBJS = $(ODIR)/main.o $(ODIR)/stockage_creux.o $(ODIR)/quadrature.o $(ODIR)/recup_node.o $(ODIR)/numerics.o
# Liste des fichiers mod
MODS = stockage_matrice.mod quadrature_P1.mod to_vtk.mod numerics.mod

# nom de l'executable
PROG = elements_finis

# regle de compilation pour le fichier executable final
$(PROG) : $(OBJS)
	$(CC) $(FFLAGS) $(OBJS) -o $(EXEDIR)/$(PROG)

# regle clean pour supprimer les fichiers binaires et executable
clean:
	rm $(ODIR)/*.o

# Regles construction des fichiers binaires
$(ODIR)/main.o: main.f90 $(ODIR)/stockage_creux.o $(ODIR)/quadrature.o $(ODIR)/recup_node.o $(ODIR)/numerics.o
	$(FC) -c $(FFLAGS) main.f90 -o $(ODIR)/main.o
$(ODIR)/numerics.o: numerics.f90 
	$(FC) -c $(FFLAGS) numerique.f90 -o $(ODIR)/numerique.o
$(ODIR)/stockage_creux.o: stockage_creux.f90 $(ODIR)/quadrature.o $(ODIR)/numerics.o
	$(FC) -c $(FFLAGS) stockage_creux.f90 -o $(ODIR)/stockage_creux.o
$(ODIR)/quadrature.o: quadrature.f90 $(ODIR)/numerics.o
	$(FC) -c $(FFLAGS) quadrature.f90 -o $(ODIR)/quadrature.o
$(ODIR)/recup_node.o: recup_node.f90 $(ODIR)/numerics.o
	$(FC) -c $(FFLAGS) recup_node.f90 -o $(ODIR)/recup_node.o