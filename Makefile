STANDART=-std=c++20
OPTIMIZATION=-O3
ANALYZER=-fanalyzer
WARNING_LVL=-Wall -Wextra -Wshadow -Wconversion -Wpedantic -Werror
BUILD_DIR=target
GCC_BUILD_DIR=${BUILD_DIR}/gcc
EXAMPLES_DIR=src/examples


examples:
	mkdir -p ${GCC_BUILD_DIR}
	g++-11 ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${EXAMPLES_DIR}/3d_fcc_Heisenberg.cpp ${SOURCES} -o ${GCC_BUILD_DIR}/3d_fcc_Heisenberg.execute
	g++-11 ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${EXAMPLES_DIR}/2d_square_Ising.cpp ${SOURCES} -o ${GCC_BUILD_DIR}/2d_square_Ising.execute

clean:
	rm -rf target
