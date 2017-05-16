# call : make "name of program"
FC=gfortran
#FFLAG1= -mcmodel=medium -ffixed-line-length-none -fbounds-check
FFLAG1= -mcmodel=medium -ffixed-line-length-none
libsac=/home/junxie/opt/sac/lib/sacio.a 
#libsac=/Users/junxie/opt/sac/lib/sacio.a 
saclib1=/home/junxie/opt/sac/lib/libsac.a
#saclib1=/Users/junxie/opt/sac/lib/libsac.a
sources=double_station.f realft.f getalpha.f gaussfilter.f four1.f getper.f distance.f envelope.f calculate_dif_az.f
objects=$(sources:.f=.o)
executable=double_station
all: $(sources) $(executable) 
.f.o:
	$(FC) $(FFLAG1) $< -c

$(executable):$(objects)
	$(FC) $(FFLAG1) $(objects) $(libsac) $(saclib1) -o $@
install:
	cp $(executable) ~/bin
uninstall:
	-rm ~/bin/getzh
clean:
	-rm $(objects) 
