#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream> 
#include <sstream>
#include <vector> 

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

        string name;
        int index;
        int nb_msat;
        int* msat_id;
        
    
    public:

        Contig() {
            this->index = -1;
            this->nb_msat = 0;
            this->name = "";

            #pragma acc enter data create(this)
            #pragma acc update device(this)
        }

        Contig(int _index, int _nb_msat, string _name) {
            this->index = _index;
            this->nb_msat = _nb_msat;
            this->name = _name;

            this->msat_id = new int[this->nb_msat];

            #pragma acc enter data create(this)
            #pragma acc update device(this)
            #pragma acc enter data create(this->msat_id[0:this->nb_msat]) 
        }

        ~Contig() {
            delete [] msat_id;
            #pragma acc exit data delete(this->msat_id[0:this->nb_msat])
            #pragma acc exit data delete(this)
        }

        int get_index() {
            return this->index;
        }

        int get_nb_msat() {
            return this->nb_msat;
        }

        string get_name() {
            return this->name;
        }

};

class Dataset {

    private:
        int nb_indiv;
        string* indiv_names;
        int nb_contigs;
        Contig** ctg;
        int nb_msat;

        void set_indiv_number(string line) {
            istringstream is(line);
            string part;
            int header = 0;
            
            while (getline(is, part, '\t')) {
                if(header > 2) {
                    this->nb_indiv++;
                }
                header++;
            }
        }

        void set_indiv_name(string line) {
            istringstream is(line);
            string part;
            int header = 0;
        
            while (getline(is, part, '\t')) {
                if(header > 2) {
                    //cout << "part : " << part << endl;
                    //cout << header << endl;
                  
                    this->indiv_names[header - 3] = part;
                }
                header++;
            }
            
            #pragma acc update device(this->indiv_names[0:this->nb_indiv])   
        }

    public:

        Dataset() {

            this->nb_indiv = 0;
            //this->indiv_names = new string[0];
            this->nb_contigs = 0;
            this->nb_msat = 0;
            
            #pragma acc enter data create(this)
            #pragma acc update device(this)
            
        }



        ~Dataset() {
        
            for (int i=0; i < this->nb_msat; ++i) {
                delete ctg[i];
            }

            delete ctg;
            #pragma acc exit data delete(this->ctg[0:this->nb_msat])
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

                    this->set_indiv_number(line);

                    indiv_names = new string[this->nb_indiv];

                    #pragma acc enter data create(this->indiv_names[0:this->nb_indiv])
                    this->set_indiv_name(line);
                }
                else {
                    this->nb_msat++;
                }

                lig++;

            }

            infile.close();

            
            if(this->nb_indiv == 0) {
                
                cout << "Error nb indiv = 0" << endl;
                exit(0);

            }           


            infile.open(file);

            //TODO need refactor

            if (infile.is_open()) {
                cout << "\nreading data from file : " << file << "\n" << endl;
            }
            else {
                cout << "\ncannot read file : " << file << endl;
                exit(0);
            }

            lig = 1;

            this->ctg = new Contig*[this->nb_msat];

            for (int i=0; i < this->nb_msat; ++i) {

                this->ctg[i] = new Contig();
                
            }

            #pragma acc enter data create(this->ctg[0:this->nb_msat])
            
            for( string line; getline( infile, line ); ) {
                
                if(lig != 1) {
                        
                }

                lig++;
            }


            cout << "reading input file is finish" << "\n" << endl;
            infile.close();    

            #pragma acc update device(this->nb_indiv)    
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

        void data_summary() {
            
            cout << this->nb_indiv << " individuals found\n" << endl;

            for (int i=0; i < this->nb_indiv; i=i+1) {
                cout << this->indiv_names[i] << " ";
            }
            cout << "\n";

            cout << this->nb_msat << " microsatellites foud in " << this->nb_contigs << endl;

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
    data->data_summary();

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
