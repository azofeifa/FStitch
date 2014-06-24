
#***************************************************************************/
#**    Author: Joey Azofeifa                                               */
#**    Date: 2/1/2014	                                                   */
#**    Description: Makefile for Segmentation.                             */
#***************************************************************************/
CXX      	= g++
CXXFLAGS 	=   -Wno-unused-variable -Wno-non-virtual-dtor -std=c++11
EXEC     	= ${PWD}/segment 	
OBJS		= ${PWD}/main.o		
SOURCE 		= ${PWD}/main.cpp 	

${EXEC}: ${OBJS} ${SOURCE}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${EXEC} 

main.o: ${PWD}/main.cpp
	${CXX} -c ${CXXFLAGS} ${PWD}/main.cpp
clean:
	rm -rf *.o 
