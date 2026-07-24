// This file is part of MARTY.
//
// MARTY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MARTY is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MARTY. If not, see <https://www.gnu.org/licenses/>.

#pragma once

#include <ostream>
#include <string>

namespace csl::libdata {

inline void replaceIntegrationMarker(std::string &text,
                                     std::string const &marker,
                                     std::string const &value)
{
    std::size_t pos = 0;
    while ((pos = text.find(marker, pos)) != std::string::npos) {
        text.replace(pos, marker.size(), value);
        pos += value.size();
    }
}

inline void print_libintegration_hdata(std::ostream &out,
                                       std::string const &libraryNamespace)
{
    std::string text = R"MTYNUM(#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "callable.h"
#include "kinematics.h"

namespace @NAMESPACE@ {

typedef std::vector<std::pair<double, double>> BoundSet;

template<class T>
struct Estimate {
    T value;
    T error;

    bool is_compatible(const Estimate<T> &other, int n_sigma=1) {
        if (this->value < other.value) {
            return this->value + n_sigma * this->error > other.value - n_sigma * other.error;
        } else {
            return this->value - n_sigma * this->error < other.value + n_sigma * other.error;
        }
    }
};

class Process {

public:

    /**
     * @brief Initializes a Process with a given squared amplitude and kinematics
     * @param f Squared amplitude function
     * @param kinematics Kinematics of the process
     */
    Process(std::unique_ptr<Callable<complex_t, param_t>> &&f, std::unique_ptr<Kinematics> &&kinematics);

    /**
     * @brief Initializes a Process with a given squared amplitude and kinematics
     * @param f_name Name of the squared amplitude function
     * @param kinematics Kinematics of the process
     */
    Process(const std::string& f_name, std::unique_ptr<Kinematics> &&kinematics);

    /**
     * @brief Semi-initializes a Process with a given kinematics and no squared amplitude function
     * @param kinematics Kinematics of the process
     */
    Process(std::unique_ptr<Kinematics> &&kinematics);

    /**
     * @brief Sets the squared amplitude of the process
     * @param f_name Name of the squared amplitude function
     */
    void set_function(const std::string& f_name);

    /**
     * @brief Overload call operator to compute the differential cross-section/decay rate of the process at a particular phase space point with possible cuts
     * @param kin_args Phase space point in relevant coordinates (see KinematicsCalculator for the relevant coordinates of each phase space)
     * @return The differential cross-section/decay rate of the process at the given phase space point
     */
    double operator()(const std::vector<double> &kin_args = {});

    /**
     * @brief Checks whether the process is fully initialized
     * @return True if the process is fully initialized, false otherwise
     */
    bool initialized();

    /**
     * @brief Yields a raw pointer to the kinematics of the process
     * @return Raw pointer to the kinematics of the process
     */
    Kinematics* get_kinematics();

private:
    /**
     * @brief Squared amplitude function
     */
    std::unique_ptr<Callable<complex_t, param_t>> m_f;

    /**
     * @brief Kinematics of the process
     * @sa Kinematics
     */
    std::unique_ptr<Kinematics> m_kinematics;

};

struct MobiusParams {
    Process* p;
    bool* has_peak;
};

class Integrator {

public:

    /**
     * @brief Partially initializes a default integrator with no process to integrate
     */
    Integrator();

    /**
     * @brief Fully initializes an integrator with a given process
     * @param process Process to be integrated
     */
    Integrator(std::unique_ptr<Process> &&process);

    /**
     * @brief Fully initializes an integrator with a default process without a squared amplitude
     * @param kinematics Kinematics of the process
     */
    Integrator(Kinematics& kinematics);

    /**
     * @brief Fully initializes an integrator with a complete process
     * @param function Name of the process' squared amplitude function
     * @param kinematics Kinematics of the process
     */
    Integrator(const std::string& function, Kinematics& kinematics);

    /**
     * @brief Sets the integrator's process
     * @param process Process to be integrated
     */
    void set_process(std::unique_ptr<Process> &&process);

    /**
     * @brief Sets the integrator's process' squared amplitude function
     * @param f_name Name of the process' squared amplitude function
     */
    void set_function(const std::string& f_name);

    /**
     * @brief Yields a pointer to the kinematics of the process
     * @return Pointer to the kinematics of the process
     */
    Kinematics* get_kinematics();

    /**
     * @brief Integrates the process using the VEGAS Monte-Carlo algorithm
     */
    void integrate();

    /**
     * @brief Integrates a process with a known sharp peak in the kinematics distribution using a Möbius transform to flatten the peaks
     * @param peak_coords Coordinates of the peak
     * @param peak_width Approximate width of the peak
     */
    void integrate_peaked(double peak_coords[], double peak_width);

    /**
     * @brief Integral of the process' with its error
     * @return An Estimate of the process' integral
     */
    Estimate<double> get_integral() const;

    /**
     * @brief Value of the process' integral
     * @return The estimated value of the process' integral
     */
    double get_integral_value() const;

    /**
     * @brief Error of the process' integral
     * @return The estimated error of the process' integral
     */
    double get_integral_error() const;

    /**
     * @brief Sets the number of points used by each iteration of VEGAS
     * @param calls Number of points
     */
    void set_calls_per_iter(size_t calls);

    /**
     * @brief Sets the maximum number of iterations of VEGAS before convergence
     * @param max_iter Maximum number of iterations
     */
    void set_max_iter(size_t max_iter);

private:
    static constexpr size_t DEFAULT_CALLS_PER_ITER = 1000;
    static constexpr size_t DEFAULT_MAX_ITER = 50;

    std::unique_ptr<Process> m_proc;
    Estimate<double> m_integral;
    bool m_converged;

    size_t m_calls_per_iter;
    size_t m_max_iter;

    /**
     * @brief Proxy function with the correct prototype for the GSL VEGAS algorithm
     * @param x Variables of integration
     * @param dim Dimension of the integration space (phase space)
     * @param par Parameters of the function (here the Process will be passed as parameter)
     * @return The differential cross-section/decay rate for the integrator's process at the given phase-space point
     */
    static double func(double x[], size_t dim, void* par);

    /**
     * @brief Proxy function after Möbius transformation of its variables with the correct prototype for the GSL VEGAS algorithm
     * @param x Variables of integration
     * @param dim Dimension of the integration space (phase space)
     * @param par Parameters of the function (here the Process will be passed as parameter)
     * @return The differential cross-section/decay rate for the integrator's process at the given phase-space point
     */
    static double mobius_func(double x[], size_t dim, void* par);

    /**
     * @brief VEGAS integration of the process' squared amplitude within the given phase space bounds
     * @param xl Lower bounds
     * @param xu Upper bounds
     * @param dim Dimension of the phase space
     * @param add If true, the integration result will be added to the current estimate of the integral instead of replacing it
     */
    void integrate_V(double *xl, double *xu, size_t dim, bool add = false);


    void integrate_MV();

    /**
     * @brief Splits the phase space along one dimension to isolate a possible peak
     * @param bounds Current sub phase spaces
     * @param x_peak Coordinate of the peak along the given dimension
     * @param width Approximate width of the peak
     * @param dim Dimension along which to split the phase space
     * @return A new set of sub-phase spaces splitted along the given dimension with the peak isolated
     */
    std::vector<BoundSet> trisect(std::vector<BoundSet>, double x_peak, double width, size_t dim);

};

} // namespace @NAMESPACE@
)MTYNUM";
    replaceIntegrationMarker(text, "@NAMESPACE@", libraryNamespace);
    out << text;
}

inline void print_libintegration_cppdata(std::ostream &out,
                                         std::string const &libraryNamespace,
                                         std::string const &groupIncludes,
                                         std::string const &functionLookup)
{
    std::string text = R"MTYNUM(#include "integration.h"

@GROUP_INCLUDES@

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace @NAMESPACE@ {


Process::Process(std::unique_ptr<Callable<complex_t, param_t>> &&f, std::unique_ptr<Kinematics> &&kinematics)
    : m_f(std::move(f)), m_kinematics(std::move(kinematics)) {}

Process::Process(const std::string& f_name, std::unique_ptr<Kinematics> &&kinematics) : m_kinematics(std::move(kinematics)) {
    set_function(f_name);
}

Process::Process(std::unique_ptr<Kinematics> &&kinematics) : m_kinematics(std::move(kinematics)) {
    m_f = nullptr;
}

void Process::set_function(const std::string &f_name) {
@FUNCTION_LOOKUP@
    throw std::out_of_range(
        "No complex-valued generated function named \"" + f_name + "\" was found.");
}

double Process::operator()(const std::vector<double>& kin_args) {
    m_kinematics->update(kin_args, true); // Update the kinematics at the given PS point
    if (m_kinematics->is_point_valid()) { // Check for possible cuts
        complex_t res = m_kinematics->get_phase_space_factor(kin_args) * (*m_f)(m_kinematics->get_params());
        if (std::abs(res.imag() / res.real()) > 1e-2) // Complex values often indicate a problem in the squared amplitude
            std::cout << "Warning: Function evaluation yielded nonzero imaginary part (" << res.imag() / res.real() << "). Check results." << std::endl;
        return res.real();
    } else {
        return 0;
    }
}

bool Process::initialized() {
    return m_kinematics && m_f;
}

Kinematics *Process::get_kinematics() {
    return m_kinematics.get();
}

Integrator::Integrator()
    : m_proc(nullptr),
      m_integral{0., 0.},
      m_converged(false),
      m_calls_per_iter(Integrator::DEFAULT_CALLS_PER_ITER),
      m_max_iter(Integrator::DEFAULT_MAX_ITER)
{}

Integrator::Integrator(std::unique_ptr<Process> &&process) : Integrator() {
    set_process(std::move(process));
}

Integrator::Integrator(Kinematics &kinematics)
: Integrator(std::make_unique<Process>(std::make_unique<Kinematics>(std::move(kinematics)))) {}

Integrator::Integrator(const std::string &function, Kinematics &kinematics)
: Integrator(std::make_unique<Process>(function, std::make_unique<Kinematics>(std::move(kinematics))))
{}

void Integrator::set_process(std::unique_ptr<Process> &&process) {
    m_proc = std::move(process);
    m_converged = false; // Reset integral
    m_integral = {0, 0};
}

void Integrator::set_function(const std::string &f_name) {
    m_proc->set_function(f_name);
    m_converged = false;
}

Kinematics *Integrator::get_kinematics() {
    return m_proc->get_kinematics();
}

void Integrator::integrate() {
    // Ensure process is fully initialized
    if (!m_proc->initialized()) {
        std::cerr << "Process to integrate is not fully initialized.\n";
        exit(123);
    }

    // Get integration bounds
    const auto bounds = m_proc->get_kinematics()->kin_limits();
    const size_t dim = m_proc->get_kinematics()->get_phase_space_dim();

    // 2-body decay phase-space is zero-dimensional, no need for integration
    if (dim == 0) {
        m_integral.value = (*m_proc)();
        m_integral.error = 0;
        m_converged = true;
        return;
    }

    std::vector<double> xl, xu;
    for (auto &&p : bounds) {
        xl.emplace_back(p.first);
        xu.emplace_back(p.second);
    }

    // Call VEGAS
    integrate_V(&xl[0], &xu[0], dim);
}

void Integrator::integrate_peaked(double peak_coords[], double peak_width) {
    // Ensure process is fully initialized
    if (!m_proc->initialized()) {
        std::cerr << "Process to integrate is not fully initialized.\n";
        exit(123);
    }

    // Get integration bounds
    const auto full_bounds = m_proc->get_kinematics()->kin_limits();
    const size_t dim = m_proc->get_kinematics()->get_phase_space_dim();

    // Split the phase space to isolate the peaks
    std::vector<std::vector<std::pair<double, double>>> bounds;
    bounds.emplace_back(full_bounds);
    for (size_t i = 0; i < dim; i++) {
        bounds = trisect(bounds, peak_coords[i], peak_width, i);
    }

    // Integrate each sub-phase space with the right algorithm
    std::vector<bool> has_peak(dim, false);
    for (BoundSet b : bounds) {
        for (size_t i = 0; i < dim; i++) {
            has_peak[i] = (b.at(i).first < peak_coords[i] && peak_coords[i] < b.at(i).second);
        }

        integrate_MV();
    }
}

Estimate<double> Integrator::get_integral() const {
    if (!m_converged) {
        std::cerr << "Integral has not been evaluated or evaluation has not converged."  << std::endl;
        exit(123);
    }
    return m_integral;
}

double Integrator::get_integral_value() const {
    if (!m_converged) {
        std::cerr << "Integral has not been evaluated or evaluation has not converged."  << std::endl;
        exit(123);
    }
    return m_integral.value;
}

double Integrator::get_integral_error() const {
    if (!m_converged) {
        std::cerr << "Integral has not been evaluated or evaluation has not converged."  << std::endl;
        exit(123);
    }
    return m_integral.error;
}

void Integrator::set_calls_per_iter(size_t calls) {
    m_calls_per_iter = calls;
}

void Integrator::set_max_iter(size_t max_iter) {
    m_max_iter = max_iter;
}

double Integrator::func(double x[], size_t dim, void *par) {
    std::vector<double> kin_p(x, x + dim);
    Process* proc = (Process*) par;
    return (*proc)(kin_p);
}

void Integrator::integrate_V(double *xl, double *xu, size_t dim, bool add) {
    double res, err;
    size_t iter {0};

    gsl_monte_function monte_func;
    monte_func.f = &func;
    monte_func.dim = dim;
    monte_func.params = m_proc.get();

    // Setup the generator
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);

    // Warm up the grid with a small run
    gsl_monte_vegas_integrate (&monte_func, xl, xu, dim, m_calls_per_iter / 10, r, s, &res, &err);
    // Iterate until convergence (chi² < 1.4)
    do {
        gsl_monte_vegas_integrate (&monte_func, xl, xu, dim, m_calls_per_iter, r, s, &res, &err);
        ++iter;
    } while (std::abs(gsl_monte_vegas_chisq(s) - 1) > 0.4 && iter < m_max_iter);
    gsl_monte_vegas_free (s);
    gsl_rng_free(r);

    if (iter == m_max_iter) {
        std::cout << "Maximum number of iterations reached. Integrator counldn't converge." << std::endl;
        m_converged = false;
    } else {
        m_converged = true;
        if (add) {
            m_integral.value += res;
            m_integral.error += err;
        } else {
            m_integral.value = res;
            m_integral.error = err;
        }
    }
}

void Integrator::integrate_MV() {}

std::vector<BoundSet> Integrator::trisect(std::vector<BoundSet> bounds, double x_peak, double width, size_t dim) {
    std::vector<BoundSet> new_bounds;
    for (auto b : bounds) {
        if (b.at(dim).first < x_peak && b.at(dim).second > x_peak) {
            new_bounds.emplace_back(b);
            new_bounds.back().at(dim) = {b.at(dim).first, x_peak - std::sqrt(width)};
            new_bounds.emplace_back(b);
            new_bounds.back().at(dim) = {x_peak - std::sqrt(width), x_peak + std::sqrt(width)};
            new_bounds.emplace_back(b);
            new_bounds.back().at(dim) = {x_peak + std::sqrt(width), b.at(dim).second};
        } else {
            new_bounds.emplace_back(b);
        }
    }
    return new_bounds;
}

} // namespace @NAMESPACE@
)MTYNUM";
    replaceIntegrationMarker(text, "@NAMESPACE@", libraryNamespace);
    replaceIntegrationMarker(text, "@GROUP_INCLUDES@", groupIncludes);
    replaceIntegrationMarker(text, "@FUNCTION_LOOKUP@", functionLookup);
    out << text;
}

} // namespace csl::libdata
