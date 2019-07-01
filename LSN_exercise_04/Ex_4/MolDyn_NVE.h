/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double estimate_pot, estimate_kin, estimate_E_tot, estimate_temp;

int run=1;      //  keeps track of how many times I have run the program.

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],x_old[m_part],y_old[m_part],z_old[m_part];
double x_new[m_part], y_new[m_part], z_new[m_part];
double v_x[m_part],v_y[m_part],v_z[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

double rescale_temperature=-666.;

// simulation
int nstep, iprint, seed;
double delta;
double sigma, epsilon, mass;
double const boltzmann = 1.38064852*1E-23;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void Conf_XYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Block_V(void);
double Block_T(void);
double Block_E_tot(void);
double Block_temp(void);
double Block_P(void);

double error(double* , double* , int );

void Previous_run(void);
void Save_for_next_run(void);

//
std::string element;
std::string state;

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
