#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream> 
#include <sstream>

/**
 * @author Falindir
 * @version 0.0.0.1
 */

using namespace std;

extern const int MAXLLINE    = 100000;
extern const int MAXLWORD    = 1000;
extern const int MAXMOTIFLG  = 10;
extern const int MAXNBALLELE = 10;

class Microsat {


};


class Contig {
    
    private:

        //string name;
        int index;
        int nb_msat;
        //int* msat_id;
        
        


    public:

        Contig() {
            this->index = -1;
            this->nb_msat = 0;

            #pragma acc enter data create(this)
            #pragma acc update device(this)
        }

        Contig(int _index, int _nb_msat) {
            this->index = _index;
            this->nb_msat = _nb_msat;

            #pragma acc enter data create(this)
            #pragma acc update device(this)
        }

        ~Contig() {
            #pragma acc exit data delete(this)
        }

        int get_index() {
            return this->index;
        }

        int get_nb_msat() {
            return this->nb_msat;
        }

};

class Dataset {

    private:
        int nb_indiv;

        int nb_contigs;

        int nb_msat;

    public:

        Dataset() {

            this->nb_indiv = 0;
            this->nb_contigs = 0;
            this->nb_msat = 0;
            
            #pragma acc enter data create(this)
            #pragma acc update device(this)
        }



        ~Dataset() {
            #pragma acc exit data delete(this)
        }

        void read_data(char* file) {
            ifstream infile (file);

            if (infile.is_open()) {
                cout << "\nreading data from file : " << file << "\n" << endl;
            }
            else {
                cout << "\ncannot read file : " << file << endl;
                exit(0);
            }

            int lig = 1;
            
            for( string line; getline( infile, line ); ) {

                
                if(lig == 1) { // read the header of file 

                    istringstream is(line);
                    string part;

                    int header = 0;

                    while (getline(is, part, '\t')) {
                        if(header > 2) { // the 3 first colone is ignored
                            cout << part << endl;
                            this->nb_indiv++;    
                        }
                        header++;
                    }
                }

                lig++;

            }
            

            cout << "reading input file is finish" << "\n" << endl;
            infile.close();       
        }

        int get_nb_indiv() {
            return this->nb_indiv;
        }

        int get_nb_contigs() {
            return this->nb_contigs;
        }
        
        int get_nb_msat() {
            return this->nb_msat;
        }
};


void usage() {

    cout << "\nusage: amp infile outfile min_cov_tot min_cov_var max_cov_tot\n" << endl;
    exit(0);
};

int tranform_param(string par) {
    
    //TODO need check
    return atoi(par.c_str());
}


int main(int argc, char *argv[]) {

    if(argc!=6) {
        usage();
    }
    
    int min_cov_ind = tranform_param(argv[3]);
    int min_cov_var = tranform_param(argv[4]);
    int max_cov_ind = tranform_param(argv[5]);

    Dataset* data = new Dataset();

    data->read_data(argv[1]);

    cout << "calculing..." << endl;

    ofstream outfile (argv[2]);

    if(outfile.is_open()) {
        cout << "\nwriting to file : " << argv[2] << "\n" << endl;

        outfile << "test 0\n"; 
    }
    else {
        cout << "cannot read file : " << argv[1] << endl;
        exit(0);
    }

    outfile.close();


    cout << "DONE\n" << endl;

}
