CXX = icl -EHsc -O2 -Qmkl -DMKL_ILP64 -Fe$@
all: run
run: percol.exe ; ./$<  && gnuplot mkpng.gp
%.exe: %.cpp ; $(CXX) $<
clean: ; -rm *.exe
