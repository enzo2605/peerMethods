all: makeall

# Compiling BLAS, CBLAS and peerMethods libraries
makeall:
	( cd external_libs/BLAS && make all)
	( cd external_libs/CBLAS && make all)
	( cd src && make all)

# Cleaning all
cleanall:
	( cd external_libs/BLAS && make clean)
	( cd external_libs/CBLAS && make cleanall)
	( cd src && make cleaner)