# makefile for ACED/GaSP

aced     = aced.o acedeval.o acedlhs.o acedoptd.o
gasp     = gasp.o gaspcv.o gaspfit.o gasppred.o gaspvis.o
crit     = crit.o critcens.o critcov.o critd.o critg.o \
        critmaxd.o critmind.o critrff.o critutil.o
database = db.o dbmanip.o dbmat.o dbmatcom.o dbmatleg.o dbscalar.o
design   = desall.o desfed.o deslhs.o desseq.o desutil.o
kriging  = krcor.o kriging.o krmatern.o krmle.o krpowexp.o krpred.o
lib      = liballoc.o libbufin.o libfile.o libin.o liblist.o libmath.o \
        libout.o libperm.o libprob.o librandn.o libreg.o \
        libsort.o libstr.o libtempl.o libvec.o
matrix   = matalloc.o matblas.o matcopy.o mateig.c matio.o matqr.o \
        matsym.o mattri.o matutil.o
minimize = min.o mincont.o minone.o minpow.o minsimp.o minxtrap.o
model    = model.o modfn.o modparse.o

# Make executables (lm is math library;
# only standard-library functions are used).

gasp: $(gasp) run.o dumcrit.o $(database) $(kriging) \
                $(lib) $(matrix) $(minimize) $(model)
	gcc $(gasp) run.o dumcrit.o $(database) $(kriging) \
		$(lib) $(matrix) $(minimize) $(model) -o gasp -lm

# Implicit rule for compiling .c to .o files.
.c.o:
	gcc -c -Wall $<
