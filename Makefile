CXX=g++-10
STANDART=-std=c++20
OPTIMIZATION=-O3
ANALYZER=-fanalyzer
WARNING_LVL=-Wall -Wextra -Wshadow -Wconversion -Wpedantic -Werror
BUILD_DIR=target
SOURCE_DIR=src
SOURCE=2dIsing.cpp
TARGET=main.execute

all:
	mkdir -p ${BUILD_DIR}
	${CXX} ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${SOURCE_DIR}/${SOURCE} -o ${BUILD_DIR}/${TARGET}

anylize:
	${CXX} ${STANDART} ${ANALYZER} ${OPTIMIZATION} ${SOURCE_DIR}/${SOURCE} -o ${BUILD_DIR}/${TARGET}