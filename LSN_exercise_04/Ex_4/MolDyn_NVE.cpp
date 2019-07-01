/*
    ***********************
    ––––– EXERCISE 4 ––––––
    ***********************
 
    Leonardo Alchieri, 886810
 
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Laboratorio di Simulazione Numerica
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
//
//
//
//
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include "MolDyn_NVE.h"
#include "random_generator/random.h"




// I use these libraries to generate a directory inline
//#include <conio.h>
//#include <dir.h>
//#include <process.h>
//#include <stdio.h>

using namespace std;

int main(int argc, char** argv){
    
    //
    //  I input the type of element I want to work on
    //
    
    
    if(argc == 1){
        element = "General";
        boltzmann = 1.;
    }
    else{
        element = argv[1];
    }
    
    if(argc < 3 ){
        state = "liquid";
    }
    else{
        state = argv[2];
    }
    
    cout << "––– The code will be run for the element: " << element << ", in the state: " << state << endl;
    
    
    
    ifstream read_run;
    read_run.open("run.dat");
    
    
    read_run >> run;
    
    char answer;
    char answer_2;
    
    //
    //  Initialize the system
    //
    if(run == 0){
        cout << "––– NEW RUN. " << endl;
        Input();
        
    }
    int patience = 0;
    
    //
    //  I check if to restart from run 0 or the previous run.
    //
    if(run != 0){
        cout << "——— The program has been run previously. Do you wish to continue from where it has been left? (y/n)" << endl;
        do{
            cin >> answer;
            if(answer == 'y'){
                cout << "** Starting from previous run." << endl;
                cout << "** Current run number: " << run << endl;
                //
                //
                /*
                 
                 *** RESCALE THE TEMPERATURE ***
                 
                 I add the possibility (as a terminal question) to rescale to match a desired temperature.
                 
                 */
                //
                //
                
                
                cout << "––– Do you want to rescale to a specific temperature? (y/n)" << endl;
                
                do{
                    cin >> answer_2;
                    if(answer_2 == 'y'){
                        cout << "* Input the temperature desired: ";
                        cin >> rescale_temperature;
                    }
                    else if(answer_2 == 'n'){
                        cout << "* The system will be rescaled according to the temperature given in the input file." << endl;
                    }
                    else{
                        cout << "--- Please answer with 'y' or 'n'. TRY AGAIN." << endl;
                    }
                }while(answer_2 != 'n' && answer_2 != 'y');
                    
                //
                //  Input function when there's a restart.
                //
                
                Previous_run();
            }
            else if(answer == 'n'){
                cout << "** Restarting from scratch." << endl;
                run=0;
                Input();
            }
            else{
                cout << "--- Please answer with 'y' or 'n'. TRY AGAIN." << endl;
            }
            patience++;
            if(patience == 50){
                cout << "--- I guess you are really stubborn. Bye." << endl;
                return 1;
            }
        }while(answer != 'y' && answer!= 'n');
    }
 
    //
    //
    //
    //
    //
    //
    int n_conf = 1;
    int blocks = 100;   //number of blocks
    int block_steps = int((double)nstep/blocks);    //number of steps in each block.
    
    double sum_T[100] = {0.};
    double sum_V[100] = {0.};
    double sum_E_tot[100] = {0.};
    double sum_temp[100] = {0.};
    double sum_P[100] = {0.};
    
    double sq_sum_T[100] = {0.};
    double sq_sum_V[100] = {0.};
    double sq_sum_E_tot[100] = {0.};
    double sq_sum_temp[100] = {0.};
    double sq_sum_P[100] = {0.};
    
    int count =0;
    
    for(int istep=1; istep <= nstep; ++istep){
        Move();           //Move particles with Verlet algorithm
        if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
        if(istep%10 == 0){
            Measure();     //Properties measurement
            Conf_XYZ(n_conf);
            //
            //  Write actual configuration in XYZ format
            //  Commented to avoid "filesystem full"!
            //
            n_conf += 1;
        }
        //
        //  Run over the different blocks
        //
        sum_T[count] += Block_T();
        sum_V[count] += Block_V();
        sum_E_tot[count] += Block_E_tot();
        sum_temp[count] += Block_temp();
        sum_P[count] += Block_P();
        
        if(istep%blocks == 0){
            sum_T[count] /= block_steps;
            sq_sum_T[count] = sum_T[count] * sum_T[count];
            
            sum_V[count] /= block_steps;
            sq_sum_V[count] = sum_V[count] * sum_V[count];
            
            sum_E_tot[count] /= block_steps;
            sq_sum_E_tot[count] = sum_E_tot[count] * sum_E_tot[count];
            
            sum_temp[count] /= block_steps;
            sq_sum_temp[count] = sum_temp[count] * sum_temp[count];
            
            sum_P[count] /= block_steps;
            sq_sum_P[count] = sum_P[count] * sum_P[count];
            
            count++;
        }
        
    }
    ConfFinal();         //Write final configuration to restart
    
    
    //
    //  *** PROGRESSIVE SUMS ***
    //
    
    ofstream ave_pot;
    ofstream ave_kin;
    ofstream ave_E_tot;
    ofstream ave_temp;
    ofstream ave_pres;
    
    ave_pot.open((element + "/" + state + "/ave_epot.out").c_str());
    ave_kin.open((element + "/" + state + "/ave_ekin.out").c_str());
    ave_E_tot.open((element + "/" + state + "/ave_etot.out").c_str());
    ave_temp.open((element + "/" + state + "/ave_temp.out").c_str());
    ave_pres.open((element + "/" + state + "/ave_pres.out").c_str());
    
    
    double progressive_sums_pot[100] = {0.};  //vector with progressive sums of the block averages
    double progressive_sums_kin[100] = {0.};
    double progressive_sums_E_tot[100] = {0.};
    double progressive_sums_temp[100] = {0.};
    double progressive_sums_pres[100] = {0.};
    
    double progressive_sq_sums_pot[100] = {0.};   //as above, but of the squared averages
    double progressive_sq_sums_kin[100] = {0.};
    double progressive_sq_sums_E_tot[100] = {0.};
    double progressive_sq_sums_temp[100] = {0.};
    double progressive_sq_sums_pres[100] = {0.};

    double progressive_errors_pot[100] = {0.};   //the errors for each progressive sum
    double progressive_errors_kin[100] = {0.};
    double progressive_errors_E_tot[100] = {0.};
    double progressive_errors_temp[100] = {0.};
    double progressive_errors_pres[100] = {0.};
    
    ave_pot << "Averages " << "Errors" << endl;
    ave_kin << "Averages " << "Errors" << endl;
    ave_E_tot << "Averages " << "Errors" << endl;
    ave_temp << "Averages " << "Errors" << endl;
    ave_pres << "Averages " << "Errors" << endl;
    for(int i=0; i<blocks; i++){
        for(int j=0; j<i+1; j++){
            progressive_sums_pot[i] += sum_V[j];
            progressive_sums_kin[i] += sum_T[j];
            progressive_sums_E_tot[i] += sum_E_tot[j];
            progressive_sums_temp[i] += sum_temp[j];
            progressive_sums_pres[i] += sum_P[j];
            
            
            progressive_sq_sums_pot[i] += sq_sum_V[j];
            progressive_sq_sums_kin[i] += sq_sum_T[j];
            progressive_sq_sums_E_tot[i] += sq_sum_E_tot[j];
            progressive_sq_sums_temp[i] += sq_sum_temp[j];
            progressive_sq_sums_pres[i] += sq_sum_P[j];
        }
        progressive_sums_pot[i] = progressive_sums_pot[i] / (i+1);        //cumulative averages
        progressive_sums_kin[i] = progressive_sums_kin[i] / (i+1);
        progressive_sums_E_tot[i] = progressive_sums_E_tot[i] / (i+1);
        progressive_sums_temp[i] = progressive_sums_temp[i] / (i+1);
        progressive_sums_pres[i] = progressive_sums_pres[i] / (i+1);
        
        progressive_sq_sums_pot[i] = progressive_sq_sums_pot[i] / (i+1);  //cumulative averages of the squares
        progressive_sq_sums_kin[i] = progressive_sq_sums_kin[i] / (i+1);
        progressive_sq_sums_E_tot[i] = progressive_sq_sums_E_tot[i] / (i+1);
        progressive_sq_sums_temp[i] = progressive_sq_sums_temp[i] / (i+1);
        progressive_sq_sums_pres[i] = progressive_sq_sums_pres[i] / (i+1);
        
        progressive_errors_pot[i] = error(progressive_sums_pot, progressive_sq_sums_pot, i);
        progressive_errors_kin[i] = error(progressive_sums_kin, progressive_sq_sums_kin, i);
        progressive_errors_E_tot[i] = error(progressive_sums_E_tot, progressive_sq_sums_E_tot, i);
        progressive_errors_temp[i] = error(progressive_sums_temp, progressive_sq_sums_temp, i);
        progressive_errors_pres[i] = error(progressive_sums_pres, progressive_sq_sums_pres, i);
        
        ave_pot << (progressive_sums_pot[i] * epsilon * boltzmann)/(double)npart << " " << (progressive_errors_pot[i] * epsilon * boltzmann)/(double)npart << endl;
        ave_kin << (progressive_sums_kin[i] * epsilon * boltzmann)/(double)npart << " " << (progressive_errors_kin[i] * epsilon * boltzmann)/(double)npart << endl;
        ave_E_tot << (progressive_sums_E_tot[i] * epsilon * boltzmann)/(double)npart << " " << (progressive_errors_E_tot[i] * epsilon * boltzmann)/(double)npart << endl;
        ave_temp << (progressive_sums_temp[i] * epsilon) << " " << (progressive_errors_temp[i] * epsilon) << endl;
        ave_pres << (progressive_sums_pres[i] * epsilon * boltzmann / (sigma*sigma*sigma)) << " " << (progressive_errors_pres[i]  * epsilon * boltzmann / (sigma*sigma*sigma)) << endl;
        
    }
    
    ave_pot.close();
    ave_pot.clear();
    
    ave_kin.close();
    ave_kin.clear();
    
    ave_E_tot.close();
    ave_E_tot.clear();
    
    ave_temp.close();
    ave_temp.clear();
    
    ave_pres.close();
    ave_pres.clear();
    //
    //
    
    run++;
    Save_for_next_run();     //save r(t) and r(t-δt) for a next run
    
    read_run.close();
    read_run.clear();   //cleans memory
    
    
    //
    //  I keep track of the number of runs saving them in external file.
    //
    ofstream write_run;
    write_run.open("run.dat");
    
    write_run << run << endl;
    
    write_run.close();
    write_run.clear();
    return 0;
}


/*
    ***************************
    ––––– IMPUT FUNCTION ––––––
    ***************************

    This function is used to load all of the initial values
    for the system - mainly the positions at the initial
    time t=0 and -δt.

*/
void Input(void){ //Prepare all stuff for the simulation
    ifstream ReadInput,ReadConf;
    //double ep, ek, pr, et, vir;

    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)¹² - (1/r)⁶]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;

    //
    //  Random number generator
    //
    //  To generate random numbers, I use an external program,
    //  which is more reliable then the srand function.
    //
    //
    
    
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("random_generator/Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
    
    ifstream input("random_generator/seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    //
    //
    //  ** end initialization random object **
    //
    //
    
    ReadInput.open(("input_" + element + "_" + state + ".dat").c_str()); //Read input

    ReadInput >> temp;

    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;
    
    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart/rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol,1.0/3.0);
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> iprint;
    
    //
    //
    //
    ReadInput >> sigma;
    ReadInput >> epsilon;
    ReadInput >> mass;
    
    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl << endl;
    ReadInput.close();
    
    //
    //  *** CHECK ***
    //
    //the length of the arrays is fixed to 108; if a bigger number of points is given, the program quits with error.
    //
    //
    if(npart > m_part){
        cerr << "––– EXIT ERROR: number of particles too big for the system. ––––" << endl;
        return;
    }
    //
    //
    //
    //Prepare array for measurements
    iv = 0; //Potential energy
    ik = 1; //Kinetic energy
    ie = 2; //Total energy
    it = 3; //Temperature
    n_props = 4; //Number of observables

    //Read initial configuration
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
    }
    ReadConf.close();

    //Prepare initial velocities
    cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){
        //
        //  For each particle, I generate a random number for velocity between -0.5 and 0.5
        //
        v_x[i] = rnd.Rannyu() - 0.5;
        v_y[i] = rnd.Rannyu() - 0.5;
        v_z[i] = rnd.Rannyu() - 0.5;
        
        sumv[0] += v_x[i];
        sumv[1] += v_y[i];
        sumv[2] += v_z[i];
    }
    //
    //  I average the velocities in all of the 3 dimensions
    //
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    //
    //
    //
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
        //
        //  I substract from each velocity the average on the specific axis
        //
        v_x[i] = v_x[i] - sumv[0];
        v_y[i] = v_y[i] - sumv[1];
        v_z[i] = v_z[i] - sumv[2];

        sumv2 += v_x[i]*v_x[i] + v_y[i]*v_y[i] + v_z[i]*v_z[i];
    }
    //
    //  I average the squares of the velocities < v² >
    //
    sumv2 /= (double)npart;

    //
    //
    //
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
    for (int i=0; i<npart; ++i){
        //
        //  I scale the velocities - in order to reflect the temperature given.
        //
         v_x[i] *= fs;
         v_y[i] *= fs;
         v_z[i] *= fs;
        
        //
        //  I load in these vectors the positions at time t-δt.
        //
         x_old[i] = Pbc(x[i] - v_x[i] * delta);
         y_old[i] = Pbc(y[i] - v_y[i] * delta);
         z_old[i] = Pbc(z[i] - v_z[i] * delta);
        
        
    }
    
    rnd.SaveSeed();
    return;
}

/*
 **********************************
 ––––– PREVIOUS RUN FUNCTION ––––––
 **********************************
 
 This function loads r(t) and r(t-δt) from a previous run. From here
 
 */
void Previous_run(void){ //Prepare all stuff for the simulation
    ifstream ReadInput,ReadConf;
    //double ep, ek, pr, et, vir;
    
    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)¹² - (1/r)⁶]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;
    
    
    
    ReadInput.open(("input_" + element + "_" + state + ".dat").c_str()); //Read input
    
    ReadInput >> temp;
    
    if(rescale_temperature == -666.)
        rescale_temperature = temp;
    
    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;
    
    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart/rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol,1.0/3.0);
    cout << "Edge of the simulation box = " << box << endl;
    
    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> iprint;
    //
    //
    //
    ReadInput >> sigma;
    ReadInput >> epsilon;
    ReadInput >> mass;
    
    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl << endl;
    ReadInput.close();
    
    //
    //  *** CHECK ***
    //
    //the length of the arrays is fixed to 108; if a bigger number of points is given, the program quits with error.
    //
    //
    if(npart > m_part){
        cerr << "––– EXIT ERROR: number of particles too big for the system. ––––" << endl;
        return;
    }
    //
    //
    //
    //Prepare array for measurements
    iv = 0; //Potential energy
    ik = 1; //Kinetic energy
    ie = 2; //Total energy
    it = 3; //Temperature
    n_props = 4; //Number of observables
    
    //Read initial configuration
    cout << "Read spatial configuration from previous run: file <run_now/now_config_run_" + to_string(run) + ".dat> " << endl << endl;
    ReadConf.open("run_now/now_config_run_" + to_string(run) + ".dat");
    for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
    }
    ReadConf.close();
    
    //
    //
    //  Load positions at time t-δt.
    //
    //
    //double sumv[3] = {0.0, 0.0, 0.0};
    ReadConf.open("run_old/old_config_run_" + to_string(run) + ".dat");
    for(int i=0; i<npart; ++i){
        //
        //  I load in these vectors the positions at time t-δt.
        //
        ReadConf >> x_old[i] >> y_old[i] >> z_old[i];
        x_old[i] = x_old[i] * box;
        y_old[i] = y_old[i] * box;
        z_old[i] = z_old[i] * box;
    }
    
    //
    //  *** RESCALING TEMPERATURE ***
    //  Rescale the velocities to meet the desired temperature.
    //
    //
    double F_x[m_part], F_y[m_part], F_z[m_part];
    //
    //  Evaluating the force on particle i.
    //
    for(int i=0; i<npart; ++i){
        F_x[i] = Force(i,0);
        F_y[i] = Force(i,1);
        F_z[i] = Force(i,2);
    }

    
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){
        //
        //  I run the Verlet algorithm to evaluate the velocities and r(t+δt)
        //
        x_new[i] = Pbc( 2.0 * x[i] - x_old[i] + F_x[i] * pow(delta,2) );
        y_new[i] = Pbc( 2.0 * y[i] - y_old[i] + F_y[i] * pow(delta,2) );
        z_new[i] = Pbc( 2.0 * z[i] - z_old[i] + F_z[i] * pow(delta,2) );
        //
        //  Calculate the velocities.
        //
        //
        v_x[i] = Pbc(x_new[i] - x_old[i])/(2.0 * delta);
        v_y[i] = Pbc(y_new[i] - y_old[i])/(2.0 * delta);
        v_z[i] = Pbc(z_new[i] - z_old[i])/(2.0 * delta);
        
        //
        //  NOTE: this method works fine, but it's not the best. The professor suggests to calculate v(t + δt/2),
        //  which works better – but it's not mandatory to change this one.
        //
        
        sumv[0] += v_x[i];
        sumv[1] += v_y[i];
        sumv[2] += v_z[i];
    }
    //
    //  I average the velocities in all of the 3 dimensions
    //
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    //
    //
    //
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
        //
        //  I substract from each velocity the average on the specific axis
        //
        v_x[i] = v_x[i] - sumv[0];
        v_y[i] = v_y[i] - sumv[1];
        v_z[i] = v_z[i] - sumv[2];
        
        sumv2 += v_x[i]*v_x[i] + v_y[i]*v_y[i] + v_z[i]*v_z[i];
    }
    //
    //  I average the squares of the velocities < v² >
    //
    sumv2 /= (double)npart;
    
    //
    //
    //
    fs = sqrt(3 * rescale_temperature / sumv2);   // fs = velocity scale factor
    for (int i=0; i<npart; ++i){
        //
        //  I scale the velocities - in order to reflect the temperature given.
        //
        v_x[i] *= fs;
        v_y[i] *= fs;
        v_z[i] *= fs;
        
        //
        //  I load in these vectors the positions at time t-δt.
        //
        //
        //  WARNING: a second order expansion doesn't work.
        //
        //x_old[i] = x_new[i] - 2 * delta * v_x[i];
        //y_old[i] = y_new[i] - 2 * delta * v_y[i];
        //z_old[i] = z_new[i] - 2 * delta * v_z[i];
        x_old[i] = Pbc(x[i] - delta * v_x[i]);
        y_old[i] = Pbc(y[i] - delta * v_y[i]);
        z_old[i] = Pbc(z[i] - delta * v_z[i]);
        
        //
        //  There was a bug – Pbc() was not present for x_old
        //
    }
    
    return;
}
/*
    **************************
    ––––– MOVE FUNCTION ––––––
    **************************

    This function moves every single particle through
    a Verlet algorithm, r(t+δt) = r(t) - r(t-δt) + δt² a(t)

    The acceleration a(t) is evaluated using Newton's Law -
    calculating the gradient of the Lennard-Jones potential —— -∇V
 
    ~~~~ See better explanation in the Force() function below. ~~~~
 
*/
void Move(void){
    //
    //  Variables for new coordinates.
    //
    double x_new, y_new, z_new;
    //
    //  Arrays for the forces acting on the i-particle.
    //
    double F_x[m_part], F_y[m_part], F_z[m_part];
    //
    //  Evaluating the force on particle i.
    //
    for(int i=0; i<npart; ++i){
    F_x[i] = Force(i,0);
    F_y[i] = Force(i,1);
    F_z[i] = Force(i,2);
    }

/*
    ************************
    ––– VERLET ALGORITHM –––
    ************************
 
    Below is applied the Verlet integration scheme.
 
*/
    for(int i=0; i<npart; ++i){
        //
        //  Calculate the new positions (in Pbc).
        //
        x_new = Pbc( 2.0 * x[i] - x_old[i] + F_x[i] * pow(delta,2) );
        y_new = Pbc( 2.0 * y[i] - y_old[i] + F_y[i] * pow(delta,2) );
        z_new = Pbc( 2.0 * z[i] - z_old[i] + F_z[i] * pow(delta,2) );
        //
        //  Calculate the velocities.
        //
        v_x[i] = Pbc(x_new - x_old[i])/(2.0 * delta);
        v_y[i] = Pbc(y_new - y_old[i])/(2.0 * delta);
        v_z[i] = Pbc(z_new - z_old[i])/(2.0 * delta);

        x_old[i] = x[i];
        y_old[i] = y[i];
        z_old[i] = z[i];

        x[i] = x_new;
        y[i] = y_new;
        z[i] = z_new;
    }
return;
}
//
//
//
/*
 ***************************
 ––––– FORCE FUNCTION ––––––
 ***************************
 
 This function gives the algorithm for evaluating the
 force that acts apon a single particle due to the
 interaction with the nearest particles.
 
 In this case, a **Lennard-Jones** potential is used.
 
 The force is obviously computed as -∇V.
 
 */
double Force(int ip, int idir){
    //
    //
    //
    //
    // ip = "the particle the force is acted upon"
    //
    double f=0.0;   //force auxiliary variable
    double dvec[3], dr;     //distances in the 3 direction and the total norm.

    for (int i=0; i<npart; ++i){
        if(i != ip){
            //
            //  Evaluation of the distance (in Pbc) between the particle "ip"
            //  and the surrounding ones.
            //
            dvec[0] = Pbc( x[ip] - x[i] );
            dvec[1] = Pbc( y[ip] - y[i] );
            dvec[2] = Pbc( z[ip] - z[i] );
            
            dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
            dr = sqrt(dr);
            //
            //  check if the distance is smaller then the "circle"
            //  around the particle "ip" where the potential acts.
            //
            if(dr < rcut){
              //
              //    The gradient of the Lennard-Jones as been
              //    already evaluated analitically.
              //
              f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
            }
        }
    }

    return f;
}
/*
 *****************************
 ––––– MEASURE FUNCTION ––––––
 *****************************
 
 In here are given, as output files, the measurements of the
 properties of the system (i.e. total kinetic energy and
 total potential energy –– needed to evaluate the conservation
 of total energy in the system).
 
 */
void Measure(){
    //int bin;
    double v=0.;
    double t=0.;
    double vij=0.;
    double dx, dy, dz, dr;
    ofstream E_pot, E_kin, E_tot, Temp, Pres;

    E_pot.open((element + "/" + state + "/output_E_pot.dat").c_str(),ios::app);
    E_kin.open((element + "/" + state +"/output_E_kin.dat").c_str(),ios::app);
    Temp.open((element + "/" + state +"/output_temp.dat").c_str(),ios::app);
    E_tot.open((element + "/" + state +"/output_E_tot.dat").c_str(),ios::app);
    Pres.open((element + "/" + state +"/output_Pres.dat").c_str(),ios::app);

    v = 0.0; //reset observables
    t = 0.0;

    //
    //  ** POTENTIAL ENERGY **
    //
    //cycle over pairs of particles
    for (int i=0; i<npart-1; ++i){
        for (int j=i+1; j<npart; ++j){
            //
            //  Evaluate distance between 2 particles
            //
            dx = Pbc( x[i] - x[j] );
            dy = Pbc( y[i] - y[j] );
            dz = Pbc( z[i] - z[j] );

            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);

            if(dr < rcut){
                //
                //  Evaluate potential (Lennard-Jones) between two particles.
                //
                vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

                //  Potential Energy
                v += vij;
            }
        }
    }

    
    //
    //  ** KINETICK ENERGY **
    //
    for (int i=0; i<npart; ++i){
        //
        //  K = ½ (v_x² + v_y^2 v_z^2)
        //
        t += 0.5 * (v_x[i]*v_x[i] + v_y[i]*v_y[i] + v_z[i]*v_z[i]);
    }
    
    //
    //  ** PRESSURE **
    //
    double estimate_P;
    double wij=0.;
    double w=0.;
    
    for (int i=0; i<npart-1; ++i){
        for (int j=i+1; j<npart; ++j){
            //
            //  Evaluate distance between 2 particles
            //
            dx = Pbc( x[i] - x[j] );
            dy = Pbc( y[i] - y[j] );
            dz = Pbc( z[i] - z[j] );
            
            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);
            
            if(dr < rcut){
                //
                //  Evaluate the virial W of a Lennard-Jones between two particles.
                //
                wij = 48.0/pow(dr,12) - 24.0/pow(dr,6);
                
                // Virial
                w += wij;
            }
        }
    }
    w /= (double)npart;
    //
    //  Use average the potential and kinetic energy of the
    //  the whole system to evaluate that of a single particle.
    //
    estimate_pot = v/(double)npart; //Potential energy
    estimate_kin = t/(double)npart; //Kinetic energy
    estimate_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    estimate_E_tot = (t+v)/(double)npart; //Total energy
    estimate_P = rho * estimate_temp + 1/(3*vol) * w;  //Pressure

    E_pot << (estimate_pot * epsilon * boltzmann)/(double)npart  << endl;
    //E_pot << estimate_pot/(double)npart  << endl;
    E_kin << (estimate_kin * epsilon * boltzmann)/(double)npart  << endl;
    //E_kin << (estimate_kin)/(double)npart  << endl;
    Temp << (estimate_temp * epsilon) << endl;
    //Temp << (estimate_temp) << endl;
    E_tot << (estimate_E_tot * epsilon * boltzmann)/(double)npart << endl;
    //E_tot << (estimate_E_tot)/(double)npart << endl;
    Pres << (estimate_P * epsilon * boltzmann / (sigma*sigma*sigma)) << endl;
    //Pres << (estimate_P) << endl;
    
    
    Pres.close();
    E_pot.close();
    E_kin.close();
    Temp.close();
    E_tot.close();

    return;
}
/*
 ***********************************
 ––––– BLOCK AVERAGE FUNCTIONS ––––––
 ***********************************
 
 This function is used to evaluate the block averages, and respective
 mean standard deviations, for the whole system.
 
 */
//
//
// I define globally these variables so to be used in more then one function.
//
//
double T_blocks=0.;       // kinetic energy
double V_blocks=0.;       // total potential
double P_blocks=0.;
double temp_block=0.;

double Block_V(){
    V_blocks=0.;
    double dx, dy, dz;
    double dr;      // distance between 2 particles
    double V_ij;    // 2 particle potential
    
    for (int i=0; i<npart-1; ++i){
        for (int j=i+1; j<npart; ++j){
            //
            //  Evaluate distance between 2 particles
            //
            dx = Pbc( x[i] - x[j] );
            dy = Pbc( y[i] - y[j] );
            dz = Pbc( z[i] - z[j] );
            
            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);
            
            if(dr < rcut){
                //
                //  Evaluate potential (Lennard-Jones) between two particles.
                //
                V_ij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
                
                //  Potential Energy
                V_blocks += V_ij;
            }
        }
    }
    
    return V_blocks/(double)npart;
}
//
//
//
double Block_T(){
    T_blocks=0.;
    for (int i=0; i<npart; ++i){
        //
        //  K = ½ (v_x² + v_y^2 v_z^2)
        //
        T_blocks += 0.5 * (v_x[i]*v_x[i] + v_y[i]*v_y[i] + v_z[i]*v_z[i]);
    }
    return T_blocks/(double)npart;
}
//
//
//
double Block_E_tot(){
    double E_tot_block;
    E_tot_block = V_blocks + T_blocks;
    return E_tot_block/(double)npart;
}

double Block_temp(){
    //double temp_blocks;
    temp_block = (2./3.) * T_blocks/(double)npart;
    return temp_block;
}
//
//
//
double Block_P(){
    double dx, dy, dz;
    double dr;      // distance between 2 particles
    double W_ij;    // 2 particle potential
    double W=0;
    for (int i=0; i<npart-1; ++i){
        for (int j=i+1; j<npart; ++j){
            //
            //  Evaluate distance between 2 particles
            //
            dx = Pbc( x[i] - x[j] );
            dy = Pbc( y[i] - y[j] );
            dz = Pbc( z[i] - z[j] );
            
            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);
            
            if(dr < rcut){
                //
                //  Evaluate potential (Lennard-Jones) between two particles.
                //
                W_ij = 48.0/pow(dr,12) - 24.0/pow(dr,6);
                
                //  Potential Energy
                W += W_ij;
            }
        }
    }
    W /= (double)npart;
    P_blocks = rho * temp_block + 1/(3*vol) * W;
    return P_blocks;
}
/*
 ******************************************
 ––––– FINAL CONFIGURATION FUNCTION ––––––
 ******************************************
 
 In this function are present the outputs to external files
 for the final configuration of the system.
 
 */
void ConfFinal(void){
    ofstream WriteConf;
    
    cout << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("config.final");

    for (int i=0; i<npart; ++i){
        WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    }
    WriteConf.close();
    return;
}
/*
 ******************************************
 –––––  i-CONFIGURATION FUNCTION ––––––
 ******************************************
 
 In this function are present the outputs to external files
 for the intermediate configurations for the system.
 
 */

void Conf_XYZ(int nconf){ //Write configuration in .xyz format
    ofstream WriteXYZ;

    
    //WriteXYZ.open("frames_run_" + to_string(run) +"/config_" + to_string(nconf) + ".xyz");
    WriteXYZ.open("frames_run_" + to_string(run) + "/config_" + to_string(nconf) + ".xyz");
    WriteXYZ << npart << endl;
    WriteXYZ << "This is only a comment!" << endl;
    for (int i=0; i<npart; ++i){
        WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
    }
    WriteXYZ.close();
    return;
}


/*
 ******************************************
 ––––– SAVE FOR NETX RUN FUNCTION ––––––
 ******************************************
 
 This function saves the positions at time r(t) and r(t-δt):
 this is done to enable a future restart.
 
 NOTE: I preferred to write a new function in order to save
 to disk the different runs indipendently.
 
 */
void Save_for_next_run(void){
    
    ofstream run_out_now;
    ofstream run_out_old;
    
    cout << "Save configuration for possible restart. " <<endl;
    //
    //  I save the different runs in separate files and folders.
    //
    run_out_now.open("run_now/now_config_run_" + to_string(run) + ".dat");
    run_out_old.open("run_old/old_config_run_" + to_string(run) + ".dat");
    
    for(int i=0; i<npart; i++){
        //
        //  **** CHECK IF IT IS CORRECT ****
        //
        run_out_now << x[i]/box << " " << y[i]/box << " " << z[i]/box << endl;
        run_out_old << x_old[i]/box << " " << y_old[i]/box << " " << z_old[i]/box << endl;
    }
    
    run_out_now.close();
    run_out_old.close();
    
    return;
}

/*
 ***************************
 ––––– PBC FUNCTION ––––––
 ***************************
 
 Algorithm for period boundary conditions woth side L=box.
 Specifically, you have "pac-man effect" in the box: if the coordinate
 is outside the box, it is returned of the other side.
 
 */
double Pbc(double r){
    return r - box * rint(r/box);
}
//
//
//
/*
 ***************************
 ––––– ERROR FUNCTION ––––––
 ***************************
 
 
 */
double error(double* avg, double* sq_avg, int n){
    
    if(n==0)
        return 0;
    else
        return sqrt((sq_avg[n]-avg[n]*avg[n])/n);
}
//
//
//
/*
 ***********************
 ––––– EXERCISE 4 ––––––
 ***********************
 
 Leonardo Alchieri, 886810
 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Laboratorio di Simulazione Numerica
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
