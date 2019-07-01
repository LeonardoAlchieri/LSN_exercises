#include "mpi.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
    MPI::Init(argc,argv);
    
    int size = MPI::COMM_WORLD.Get_size(); int rank = MPI::COMM_WORLD.Get_rank();
    
    cout<<" Sono il nodo "<<rank<<" dei "<<size<<" che hai utilizzato!"<<endl;
    
    MPI::Finalize();
    return 0;
    
}
