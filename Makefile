


all: analyse_msat_polym.cc	
	pgc++ -acc -Minfo -std=c++11 -o GPU.out analyse_msat_polym.cc
	
	./GPU.out Caenorhabditis_brenneri.res result_GPU.csv 10 2 10000


