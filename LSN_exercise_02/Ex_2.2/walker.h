#ifndef __walker__
#define __walker__

class walker {
    
private:
    double x, y, z;
    double step;
protected:
    
public:
    // constructors
    walker(double step_lenght, double x_0, double y_0, double z_0);
    walker();
    walker(double step_length);
    // destructor
    ~walker();
    // methods
    double get_step();
    
    void movement_discrete(int direction);
    void movement_continuum(double theta, double phi);
    
    double get_distance();
    
    double get_x();
    double get_y();
    double get_z();
};

#endif // __walker__

