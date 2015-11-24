#***************************************************************************/
#**    Author: Joey Azofeifa                                               */
#**    Date:   1/2/2015	                                                   */
#**    Description: Makefile for FastReadStitcher                          */
#***************************************************************************/

CXX             = g++
CXXFLAGS        = -Wno-unused-variable -Wno-non-virtual-dtor -std=c++11  -fopenmp -Wno-write-strings
EXEC            = ${PWD}/FStitch
ARCH  			= getconf LONG_BIT

CPP_FLAGS_32  	= -D32_BIT
CPP_FLAGS_64  	= -D64_BIT

GCCVERSION 		= $(shell gcc --version | grep ^gcc | sed 's/^.* //g')
FStitch: main.o BaumWelch.o NewtonsMethod.o read_in_parameters.o read.o main_train.o grabTrainingExamples.o interval_tree.o \
		 split.o write.o main_segment.o viterbi.o validate.o
	@${CXX} ${CXXFLAGS} ${PWD}/main.o  ${PWD}/BaumWelch.o \
		${PWD}/NewtonsMethod.o ${PWD}/read_in_parameters.o ${PWD}/read.o \
		${PWD}/interval_tree.o ${PWD}/grabTrainingExamples.o ${PWD}/main_train.o \
		${PWD}/split.o ${PWD}/write.o ${PWD}/main_segment.o ${PWD}/viterbi.o ${PWD}/validate.o  -o ${EXEC}
	@echo "========================================="
	@echo Sucessfully Compiled
		
main.o: 
		@echo "========================="
		
		@${CXX} -c ${CXXFLAGS} ${PWD}/main.cpp 
		@echo Compiled main
read_in_parameters.o: 
		@${CXX} -c ${CXXFLAGS} ${PWD}/read_in_parameters.cpp 
		@echo Compiled read_in_parameters
read.o: 		
		@@${CXX} -c ${CXXFLAGS} ${PWD}/read.cpp 
		@echo Compiled read
main_train.o: 		
		@${CXX} -c ${CXXFLAGS} ${PWD}/main_train.cpp 
		@echo Compiled main_train
split.o: 
		@@${CXX} -c ${CXXFLAGS} ${PWD}/split.cpp
		@echo Compiled split
grabTrainingExamples.o: 
		@${CXX} -c ${CXXFLAGS} ${PWD}/grabTrainingExamples.cpp 
		@echo Compiled grabTrainingExamples
interval_tree.o: 
		@${CXX} -c ${CXXFLAGS} ${PWD}/interval_tree.cpp 
		@echo Compiled interval_tree
NewtonsMethod.o:
		@${CXX} -c ${CXXFLAGS} ${PWD}/NewtonsMethod.cpp 
		@echo Compiled NewtonsMethod
BaumWelch.o: 
		@${CXX} -c ${CXXFLAGS} ${PWD}/BaumWelch.cpp 
		@echo Compiled BaumWelch
main_segment.o:
		@${CXX} -c ${CXXFLAGS} ${PWD}/main_segment.cpp 
		@echo Compiled main_segment
viterbi.o:
		@${CXX} -c ${CXXFLAGS} ${PWD}/viterbi.cpp 
		@echo Compiled viterbi
validate.o:
		@${CXX} -c ${CXXFLAGS} ${PWD}/validate.cpp 
		@echo Compiled validate
		
write.o:
	@${CXX} -c ${CXXFLAGS} ${PWD}/write.cpp 
	@echo Compiled write
clean:
	@rm -f *.o


		

