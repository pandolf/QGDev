
CFLAGS=`root-config --libs --cflags`

.PHONY:all
all: QGLikelihoodCalculator.o AddLikelihoodFriend

QGLikelihoodCalculator.o: QGLikelihoodCalculator.C QGLikelihoodCalculator.h ReadTxt.C ReadParameters.C PtBins.h
	g++ $(CFLAGS) -c QGLikelihoodCalculator.C
AddLikelihoodFriend: AddLikelihoodFriend.C QGLikelihoodCalculator.o
	g++ -D STANDALONE $(CFLAGS) -o AddLikelihoodFriend AddLikelihoodFriend.C
