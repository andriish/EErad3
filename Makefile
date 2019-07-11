NAME1    = eerad3
NAME2    = eerad3_combine
NAME3    = eerad3_dist
NAME4    = eecanalytic

all: eerad3 eerad3_combine eerad3_dist eecanalytic
SOURCEDIR = ./src
OBJDIR = ./obj

VPATH = $(SOURCEDIR)


FFILES1   = eerad3.f histo.f ecuts.f eerad3lib.f phaseee.f sig.f aversub0.f aversub1.f aversub2.f virt.f brem.f 3jme.f tdhpl.f hplog.f
FFILES2   = eerad3_combine.f
FFILES3   = eerad3_dist.f
FFILES4   = src/eecanalytic.f 

#for gfortran compiler
FC        = gfortran
#FFLAGS    = -g -fno-automatic -O -finit-integer=0  -Wall -fcheck=all -g -fbacktrace  -finit-real=zero -ffpe-trap=invalid,zero,overflow,underflow
#FFLAGS    = -g -fno-automatic -O -finit-integer=0    -finit-real=zero -ffpe-trap=invalid,zero,overflow,underflow 
FFLAGS    = -fno-automatic -O2 -finit-integer=0    -finit-real=zero 
#for ifort compiler
#FC        = ifort
#FFLAGS    = -save -O4
LDFLAGS =  -static

OBJFILES1 = $(addprefix $(OBJDIR)/,$(patsubst %.f,%.o,$(FFILES1)))
OBJFILES2 = $(addprefix $(OBJDIR)/,$(patsubst %.f,%.o,$(FFILES2)))
OBJFILES3 = $(addprefix $(OBJDIR)/,$(patsubst %.f,%.o,$(FFILES3)))

$(OBJDIR)/%.o:	%.f
	$(FC) $(FFLAGS) -c $< -o $@

$(NAME1): $(OBJFILES1)
	$(FC) $(LDFLAGS) $(FFLAGS) -o $@ $(OBJFILES1)
$(NAME2): $(OBJFILES2)
	$(FC) $(LDFLAGS) $(FFLAGS) -o $@ $(OBJFILES2)
$(NAME3): $(OBJFILES3)
	$(FC) $(LDFLAGS) $(FFLAGS) -o $@ $(OBJFILES3) 

# We avoid -fcheck=all
FFLAGS4    = -g -fno-automatic -O -finit-integer=0    -finit-real=zero -ffpe-trap=invalid,zero,overflow,underflow  
$(NAME4): $(FFILES4)
	$(FC) $(LDFLAGS) $(FFLAGS4) $(FFILES4) -o $@ 

clean:
	rm -f obj/* eerad3 eerad3_dist eerad3_combine eecanalytic


