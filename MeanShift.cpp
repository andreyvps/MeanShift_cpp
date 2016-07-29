#include <stdio.h>
#include <math.h>
#include "MeanShift.h"

using namespace std;

double euclidean_distance(const vector<double> &point_a, const vector<double> &point_b){
    double total = 0;
    for(int i=0; i<point_a.size(); i++){
        total += (point_a[i] - point_b[i]) * (point_a[i] - point_b[i]);
    }
    return sqrt(total);
}

double gaussian_kernel(double distance, double kernel_bandwidth){
    double temp =  exp(-(distance*distance) / (kernel_bandwidth));
    return temp;
}

void MeanShift::set_kernel( double (*_kernel_func)(double,double) ) {
    if(!_kernel_func){
        kernel_func = gaussian_kernel;
    } else {
        kernel_func = _kernel_func;    
    }
}

vector<Cluster> MeanShift::cluster(
    const vector<vector<double> > & points,
    const double & kernel_bandwidth,
    const vector<vector<double> > & probes)
{
    vector<vector<double>> mutable_probes = probes;
    meanshift(points, kernel_bandwidth, mutable_probes);
    vector<vector<double>> modes = find_modes(mutable_probes);
    return cluster(points, modes);
}

void MeanShift::meanshift(
    const std::vector<std::vector<double> > & points,
    const double & kernel_bandwidth,
    std::vector<std::vector<double> > & probes) 
{
    vector<bool> stop_moving(probes.size(), false);
    double max_shift_distance;
    do {
        max_shift_distance = shift_probes(points, kernel_bandwidth, probes, stop_moving);
        printf("max_shift_distance: %f\n", max_shift_distance);
    } while (max_shift_distance > epsilon);
}

// Shift the probes to their local modes
double MeanShift::shift_probes(
    const std::vector<std::vector<double> > & points,
    const double & kernel_bandwidth,
    std::vector<std::vector<double> > & probes,
    vector<bool> & stop_moving)
{
    double max_shift_distance = 0;
    for (int i = 0; i<probes.size(); i++){
        if (!stop_moving[i]) {
            double shift_distance = shift_probe(points, kernel_bandwidth, probes[i]);
            if (shift_distance > max_shift_distance){
                max_shift_distance = shift_distance;
            }
            if (shift_distance <= epsilon) {
                stop_moving[i] = true;
            }
        }
    }
    return max_shift_distance;
}

double MeanShift::shift_probe(
    const std::vector<std::vector<double> > & points,
    const double & kernel_bandwidth,
    std::vector<double> & probe)
{
    vector<double> shifted_point(probe.size(), 0);
    double total_weight = 0;
    for(int i=0; i<points.size(); i++){
        vector<double> temp_point = points[i];
        double distance = euclidean_distance(probe, temp_point);
        double weight = kernel_func(distance, kernel_bandwidth);
        for(int j=0; j<shifted_point.size(); j++){
            shifted_point[j] += temp_point[j] * weight;
        }
        total_weight += weight;
    }

    for(int i=0; i<shifted_point.size(); i++){
        shifted_point[i] /= total_weight;
    }

    double distance = euclidean_distance(probe, shifted_point);
    for (int i = 0; i < shifted_point.size(); i++) {
        probe[i] = shifted_point[i];
    }

    return distance;
}

std::vector<vector<double> > MeanShift::find_modes(
    const std::vector<std::vector<double> > & converged_probes)
{
    std::vector<vector<double> > modes;

    for (int i = 0; i < converged_probes.size(); i++) {
        vector<double> probe = converged_probes[i];

        int m = 0;
        for (; m < modes.size(); m++) {
            if (euclidean_distance(modes[m], probe) < 2 * epsilon) {
                break;
            }
        }

        if (m == modes.size()) {
            modes.push_back(probe);
        }
    }

    return modes;
}

vector<Cluster> MeanShift::cluster(
    const vector<vector<double> > & points,
    const vector<vector<double> > & modes)
{
    vector<Cluster> clusters(modes.size());
    for (int m = 0; m < modes.size(); m++) {
        clusters[m].mode = modes[m];
    }

    for (int i = 0; i < points.size(); i++) {

        int closest_mode_index = 0;
        double distance_to_closest_mode = euclidean_distance(points[i], modes[0]);
        for (int c = 1; c < clusters.size(); c++) {
            double distance = euclidean_distance(points[i], modes[c]);
            if (distance < distance_to_closest_mode) {
                distance_to_closest_mode = distance;
                closest_mode_index = c;
            }
        }

        clusters[closest_mode_index].support.push_back(points[i]);
    }

    return clusters;
}