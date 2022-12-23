STANDART=-std=c++20
OPTIMIZATION=-O3
ANALYZER=-fanalyzer
WARNING_LVL=-Wall -Wextra -Wshadow -Wconversion -Wpedantic -Werror
BUILD_DIR=target
GCC_BUILD_DIR=${BUILD_DIR}/gcc
CLANG_BUILD_DIR=${BUILD_DIR}/clang
EXAMPLES_DIR=src/examples


examples:
	mkdir -p ${CLANG_BUILD_DIR}
	clang++ ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${EXAMPLES_DIR}/3d_fcc_Heisenberg.cpp ${SOURCES} -o ${CLANG_BUILD_DIR}/3d_fcc_Heisenberg.execute
	clang++ ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${EXAMPLES_DIR}/2d_square_Ising.cpp ${SOURCES} -o ${CLANG_BUILD_DIR}/2d_square_Ising.execute
	clang++ ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${EXAMPLES_DIR}/3d_fcc_Heisenberg_Multilayer.cpp ${SOURCES} -o ${CLANG_BUILD_DIR}/3d_fcc_Heisenberg_Multilayer.execute
	clang++ ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${EXAMPLES_DIR}/3d_fcc_Heisenberg_Multilayer_Current.cpp ${SOURCES} -o ${CLANG_BUILD_DIR}/3d_fcc_Heisenberg_Multilayer_Current.execute

	mkdir -p ${GCC_BUILD_DIR}
	g++-10 ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${EXAMPLES_DIR}/3d_fcc_Heisenberg.cpp ${SOURCES} -o ${GCC_BUILD_DIR}/3d_fcc_Heisenberg.execute
	g++-10 ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${EXAMPLES_DIR}/2d_square_Ising.cpp ${SOURCES} -o ${GCC_BUILD_DIR}/2d_square_Ising.execute
	g++-10 ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${EXAMPLES_DIR}/3d_fcc_Heisenberg_Multilayer.cpp ${SOURCES} -o ${GCC_BUILD_DIR}/3d_fcc_Heisenberg_Multilayer.execute
	g++-10 ${STANDART} ${OPTIMIZATION} ${WARNING_LVL} ${EXAMPLES_DIR}/3d_fcc_Heisenberg_Multilayer_Current.cpp ${SOURCES} -o ${GCC_BUILD_DIR}/3d_fcc_Heisenberg_Multilayer_Current.execute

debug:
	mkdir -p ${GCC_BUILD_DIR}
	g++-10 ${STANDART} ${OPTIMIZATION} -g3 ${WARNING_LVL} ${EXAMPLES_DIR}/3d_fcc_Heisenberg.cpp ${SOURCES} -o ${GCC_BUILD_DIR}/3d_fcc_Heisenberg.execute
	g++-10 ${STANDART} ${OPTIMIZATION} -g3 ${WARNING_LVL} ${EXAMPLES_DIR}/2d_square_Ising.cpp ${SOURCES} -o ${GCC_BUILD_DIR}/2d_square_Ising.execute
	g++-10 ${STANDART} ${OPTIMIZATION} -g3 ${WARNING_LVL} ${EXAMPLES_DIR}/3d_fcc_Heisenberg_Multilayer.cpp ${SOURCES} -o ${GCC_BUILD_DIR}/3d_fcc_Heisenberg_Multilayer.execute
	g++-10 ${STANDART} ${OPTIMIZATION} -g3 ${WARNING_LVL} ${EXAMPLES_DIR}/3d_fcc_Heisenberg_Multilayer_Current.cpp ${SOURCES} -o ${GCC_BUILD_DIR}/3d_fcc_Heisenberg_Multilayer_Current.execute
	
clean:
	rm -rf target
