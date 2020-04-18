#
# makefile                         
# 
     fc = gfortran
     fflags = $(myflgs)

     libg= -lm -llapack -lblas 
     
 
.SUFFIXES: .o .f .f90

.f90.o:
	$(fc) $(fflags) -c $*.f90
.f.o:
	$(fc) $(fflags) -c $*.f

obj  = FMSAVE.o FMZM90.o FM.o lanczos.o bounds_MP.o

src  = FMSAVE.f90 FMZM90.f90 FM.f90 lanczos.f90 bounds_MP.f90

lower: LowerBound.f90   $(obj)  $(src)
	$(fc) $(fflags)  -o $(HOME)/bin/lower $(obj) LowerBound.f90 $(libg)

clean:
	/bin/rm -f  *.o *.bak core *.mod *~ *.pdf *.ps

realclean:
	/bin/rm -f  *.o *.bak core *.mod *~ *.pdf *.ps $(HOME)/bin/lower*
