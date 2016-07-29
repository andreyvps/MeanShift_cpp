#pragma once 

#include <vector>

struct Cluster {
    std::vector<double> mode;
    std::vector<std::vector<double>> support;
};

class MeanShift {
private:
    double epsilon;
public:
    MeanShift() 
    : epsilon(0.0000001) 
    { set_kernel(NULL); }
    
    MeanShift(double _epsilon) 
    : epsilon(_epsilon) 
    { set_kernel(NULL); }
    
    MeanShift(double _epsilon, double (*_kernel_func)(double,double))
    : epsilon(_epsilon)
    { set_kernel(kernel_func); }
    
    // Find clusters using all points as probes
    std::vector<Cluster> cluster(
        const std::vector<std::vector<double> > & points, 
        const double & kernel_bandwidth,
        const std::vector<std::vector<double> > & probes);
    
    // Shift the probes to the modes of the distributions defined by the points 
    // and the bandwidth until they converge.
    void meanshift(
        const std::vector<std::vector<double> > & points, 
        const double & kernel_bandwidth,
        std::vector<std::vector<double> > & probes);
        
    // Shift the probes to their local modes
    double shift_probes(
        const std::vector<std::vector<double> > & points, 
        const double & kernel_bandwidth,
        std::vector<std::vector<double> > & probes,
        std::vector<bool> & stop_moving);
    
    // Shift just one probe over the surface defined by the points and bandwidth
    double shift_probe(
        const std::vector<std::vector<double> > & points, 
        const double & kernel_bandwidth,
        std::vector<double> & probe);
    
    // Probes converge to a set of modes.
    std::vector<std::vector<double> > find_modes(
        const std::vector<std::vector<double> > & converged_probes);
    
    // Create clusters by assigning points to the nearest mode.
    std::vector<Cluster> cluster(
        const std::vector<std::vector<double> > & points, 
        const std::vector<std::vector<double> > & modes);

private:
    double (*kernel_func)(double,double);
    void set_kernel(double (*_kernel_func)(double,double));
};
