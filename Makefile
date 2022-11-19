CXX=g++-10
STANDART=-std=c++20
OPTIMIZATION=-O3
ANALYZER=-fanalyzer
WARNING_LVL=-Wall -Wextra -Wshadow -Wconversion -Wpedantic -Werror
BUILD_DIR=target
SOURCE_DIR=src
TARGET=main.execute

all:
	mkdir -p ${BUILD_DIR}
	${CXX} ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${SOURCE_DIR}/3d_fcc_Heisenberg.cpp -o ${BUILD_DIR}/${TARGET}

2d_square_Ising:
	mkdir -p ${BUILD_DIR}
	${CXX} ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${SOURCE_DIR}/2d_square_Ising.cpp -o ${BUILD_DIR}/${TARGET}

anylize:
	${CXX} ${STANDART} ${ANALYZER} ${OPTIMIZATION} ${SOURCE_DIR}/${SOURCE} -o ${BUILD_DIR}/${TARGET}