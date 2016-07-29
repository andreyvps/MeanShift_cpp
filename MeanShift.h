#pragma once 

#include <vector>


namespace meanshift {

    typedef std::vector<double> Point;

    double euclidean_distance(const Point &point_a, const Point &point_b);
    double gaussian_kernel(double distance, double kernel_bandwidth);

    class MeanShift {
    private:
        double epsilon;
    public:
        MeanShift()
            : epsilon(0.0000001)
        {
            set_kernel(NULL);
        }

        MeanShift(double _epsilon)
            : epsilon(_epsilon)
        {
            set_kernel(NULL);
        }

        MeanShift(double _epsilon, double(*_kernel_func)(double, double))
            : epsilon(_epsilon)
        {
            set_kernel(kernel_func);
        }

        // Use meanshift on the probes and return the modes that they converge to
        std::vector<Point> find_modes(
            const std::vector<Point> & points,
            const double & kernel_bandwidth,
            std::vector<Point> & probes);

        // Shift the probes to the modes of the distributions defined by the points 
        // and the bandwidth until they converge.
        void meanshift(
            const std::vector<Point> & points,
            const double & kernel_bandwidth,
            std::vector<Point> & probes);

        // Shift the probes to their local modes
        double shift_probes(
            const std::vector<Point> & points,
            const double & kernel_bandwidth,
            std::vector<Point> & probes,
            std::vector<bool> & stop_moving);

        // Shift just one probe over the surface defined by the points and bandwidth
        double shift_probe(
            const std::vector<Point> & points,
            const double & kernel_bandwidth,
            std::vector<double> & probe);

        // Probes converge to a set of modes.
        std::vector<Point> find_modes(
            const std::vector<Point> & converged_probes);

    private:
        double(*kernel_func)(double, double);
        void set_kernel(double(*_kernel_func)(double, double));
    };

}