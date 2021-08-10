//! \file openmc_driver.h
//! Driver to initialize and run OpenMC in stages
#ifndef ENRICO_OPENMC_DRIVER_H
#define ENRICO_OPENMC_DRIVER_H

#include "enrico/cell_instance.h"
#include "enrico/geom.h"
#include "enrico/neutronics_driver.h"

#include "openmc/cell.h"
#include "openmc/tallies/filter_cell_instance.h"
#include "openmc/tallies/tally.h"
#include <gsl/gsl>
#include <mpi.h>

#include <vector>

namespace enrico {
class BoronDriverOpenmc : public BoronDriver {
public:
  explicit BoronDriverOpenmc(MPI_Comm comm);

  ~BoronDriverOpenmc();

  //! Gets the boron concentration and water density of the most current Openmc run and
  //! \return Boron concentrationi in [ppm]
  double get_boron_ppm();

  //! Gets the water density of the borated water
  //! \return Water density in [g/cm^3]
  double get_H2O_density() const override;

  //! Sets the current and previous keff in the boron driver class
  //! \param keff is the current k-effective after the Openmc run
  //! \param keffprev is the previous k-effective
  void set_k_effective(double keff, double keffprev);

  void set_ppm(double ppm, double ppm_prev);

  //! Estimates the boron concentration in ppm to find criticality condition
  //! \param step is used to determine if an initial slope is needed
  //! \return Boron concentration in [ppm]
  double solveppm(int step);

private:
  double k_eff_;
  double k_eff_prev;
  double ppm_;
  double ppm_prev_;
  double H2O_dens_;
};

//! Driver to initialize and run OpenMC in stages
class OpenmcDriver : public NeutronicsDriver {
public:
  //! One-time initalization of OpenMC and member variables
  //! \param comm An existing MPI communicator used to inialize OpenMC
  explicit OpenmcDriver(MPI_Comm comm);

  //! One-time finalization of OpenMC
  ~OpenmcDriver();

  //////////////////////////////////////////////////////////////////////////////
  // NeutronicsDriver interface

  //! Find cells corresponding to a vector of positions
  //! \param positions (x,y,z) coordinates to search for
  //! \return Handles to cells
  std::vector<CellHandle> find(const std::vector<Position>& position) override;

  //! Set the boron concentration of the material in a cell
  //! \param ppm Boron Concentration in [ppm]
  void set_boron_ppm(double ppm, double H2Odens);

  //! Set the density of the material in a cell
  //! \param cell Handle to a cell
  //! \param rho Density in [g/cm^3]
  void set_density(CellHandle cell, double rho) const override;

  //! Set the temperature of a cell
  //! \param cell Handle to a cell
  //! \param T Temperature in [K]
  void set_temperature(CellHandle cell, double T) const override;

  //! Get the density of a cell
  //! \param cell Handle to a cell
  //! \return Cell density in [g/cm^3]
  double get_density(CellHandle cell) const override;

  //! Get the temperature of a cell
  //! \param cell Handle to a cell
  //! \return Temperature in [K]
  double get_temperature(CellHandle cell) const override;

  //! Get the volume of a cell
  //! \param cell Handle to a cell
  //! \return Volume in [cm^3]
  double get_volume(CellHandle cell) const override;

  //! Detemrine whether a cell contains fissionable nuclides
  //! \param cell Handle to a cell
  //! \return Whether the cell contains fissionable nuclides
  bool is_fissionable(CellHandle cell) const override;

  std::size_t n_cells() const override { return cells_.size(); }

  //! Create energy production tallies
  void create_tallies() override;

  //! Determine number of cells participating in coupling
  //! \return Number of cells
  xt::xtensor<double, 1> heat_source(double power) const final;

  std::string cell_label(CellHandle cell) const;

  gsl::index cell_index(CellHandle cell) const override;

  //////////////////////////////////////////////////////////////////////////////
  // Driver interface

  //! Initialization required in each Picard iteration
  void init_step() final;

  //! Runs OpenMC for one Picard iteration
  void solve_step() final;

  double k_effective() const override;

  //! Writes OpenMC output for given timestep and iteration
  //! \param timestep timestep index
  //! \param iteration iteration index
  void write_step(int timestep, int iteration) final;

  //! Finalization required in each Picard iteration
  void finalize_step() final;

private:
  // Data members
  openmc::Tally* tally_;                     //!< Fission energy deposition tally
  openmc::CellInstanceFilter* filter_;       //!< Cell instance filter
  std::map<CellHandle, CellInstance> cells_; //!< Array of cell instances
  int n_fissionable_cells_;                  //!< Number of fissionable cells in model
};

} // namespace enrico

#endif // ENRICO_OPENMC_DRIVER_H
