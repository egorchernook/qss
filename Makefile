STANDART=-std=c++20
OPTIMIZATION=-O3
ANALYZER=-fanalyzer
WARNING_LVL=-Wall -Wextra -Wshadow -Wconversion -Wpedantic -Werror
BUILD_DIR=target
SOURCE_DIR=src/examples

examples:
	mkdir -p ${BUILD_DIR}/gcc
	mkdir -p ${BUILD_DIR}/clang
	g++-10 ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${SOURCE_DIR}/3d_fcc_Heisenberg.cpp -o ${BUILD_DIR}/gcc/3d_fcc_Heisenberg.execute
	g++-10 ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${SOURCE_DIR}/2d_square_Ising.cpp -o ${BUILD_DIR}/gcc/2d_square_Ising.execute

	# clang ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${SOURCE_DIR}/3d_fcc_Heisenberg.cpp -o ${BUILD_DIR}/clang/3d_fcc_Heisenberg.execute
	# clang ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${SOURCE_DIR}/2d_square_Ising.cpp -o ${BUILD_DIR}/clang/2d_square_Ising.execute

clean:
	rm -rf target