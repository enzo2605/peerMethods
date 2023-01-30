all: makeall

# Compiling peerMethods libraries
makeall:
	( cd src && make all)

# Cleaning all
cleanall:
	( cd src && make cleanall)