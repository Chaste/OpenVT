/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef POTTSPARABOLICPDE_HPP_
#define POTTSPARABOLICPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "PottsBasedCellPopulation.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractLinearParabolicPde.hpp"

/**
 * A parabolic PDE to be solved numerically using the finite element method, for
 * coupling to a cell-based simulation.
 *
 * The PDE takes the form
 *
 * c*du/dt = Grad.(D*Grad(u)) + (a*u + b) * rho(x),
 *
 * where the scalars c, D, a and b are specified by the members mDuDtCoefficient,
 * mDiffusionCoefficient, mLinearSourceCoefficient and mConstantSourceCoefficient, 
 * respectively. Their values must be set in the constructor.
 *
 * The function rho(x) denotes the local density of non-apoptotic cells. This
 * quantity is computed for each element of a 'coarse' finite element mesh that is
 * passed to the method SetupSourceTerms() and stored in the member mCellDensityOnCoarseElements.
 * For a point x, rho(x) is defined to be the number of non-apoptotic cells whose
 * centres lie in each finite element containing that point, scaled by the area of
 * that element.
 */
template<unsigned DIM>
class PottsParabolicPde : public AbstractLinearParabolicPde<DIM,DIM>
{
    friend class TestCellBasedParabolicPdes;

private:

    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the PDE and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<AbstractLinearParabolicPde<DIM, DIM> >(*this);
       archive & mConstantCellSourceCoefficient;
       archive & mLinearCellSourceCoefficient;
       archive & mConstantSourceCoefficient;
       archive & mLinearSourceCoefficient;
       archive & mDiffusionCoefficient;
       archive & mDuDtCoefficient;
       archive & mCellDensityOnCoarseElements;
    }

protected:

    /** The cell population member. */
    PottsBasedCellPopulation<DIM>& mrCellPopulation;

    /** Coefficient of constant density dependent source term. */
    double mConstantCellSourceCoefficient;

    /** Coefficient of linear density dependent source term. */
    double mLinearCellSourceCoefficient;

    /** Coefficient of constant source term. */
    double mConstantSourceCoefficient;

    /** Coefficient of linear source term. */
    double mLinearSourceCoefficient;

    /** Diffusion coefficient. */
    double mDiffusionCoefficient;

    /** Coefficient of rate of change term.  */
    double mDuDtCoefficient;

    /** Vector of averaged cell densities on elements of the coarse mesh. */
    std::vector<double> mCellDensityOnCoarseElements;

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation reference to the potts based cell population
     * @param constantCellSourceCoefficient the constant density dependent source term coefficient (defaults to 0.0)
     * @param linearCellSourceCoefficient the linear density dependent source term coefficient (defaults to 0.0)
     * @param constantSourceCoefficient the constant source term coefficient (defaults to 0.0)
     * @param linearSourceCoefficient the linear source term coefficient (defaults to 0.0)
     * @param diffusionCoefficient the rate of diffusion (defaults to 1.0)
     * @param duDtCoefficient rate of reaction (defaults to 1.0)
     */
    PottsParabolicPde(PottsBasedCellPopulation<DIM>& rCellPopulation,
                               double constantCellSourceCoefficient=0.0, 
                               double linearCellSourceCoefficient=0.0,
                               double constantSourceCoefficient=0.0, 
                               double linearSourceCoefficient=0.0,
                               double diffusionCoefficient=1.0,
                               double duDtCoefficient=1.0);

    /**
     * @return const reference to the cell population (used in archiving).
     */
    const PottsBasedCellPopulation<DIM>& rGetCellPopulation() const;

    /**
     * @return mConstantCellSourceCoefficient
     */
    double GetConstantCellCoefficient() const;

    /**
     * @return mLinearCellSourceCoefficient
     */
    double GetLinearCellCoefficient() const;

    /**
     * @return mConstantSourceCoefficient
     */
    double GetConstantCoefficient() const;

    /**
     * @return mLinearSourceCoefficient
     */
    double GetLinearCoefficient() const;
        
    /**
     * @return mDiffusionCoefficient
     */
    double GetDiffusionCoefficient() const;

    /**
     * @return mDuDtCoefficient
     */
    double GetDuDtCoefficient() const;

    /**
     * Set up the source terms.
     *
     * \todo this is identical to the one in PottsEllipticPde so refactor.
     *
     * @param rCoarseMesh reference to the coarse mesh
     * @param pCellPdeElementMap optional pointer to the map from cells to coarse elements
     */
    void virtual SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap=nullptr);

    /**
     * Overridden ComputeDuDtCoefficientFunction() method.
     *
     * @return the function c(x) in "c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)"
     *
     * @param rX the point in space at which the function c is computed
     */
    virtual double ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& rX);

    /**
     * Overridden ComputeSourceTerm() method.
     *
     * @return computed source term.
     *
     * @param rX the point in space at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the point
     * @param pElement the mesh element that x is contained in (optional; defaults to NULL).
     */
    virtual double ComputeSourceTerm(const ChastePoint<DIM>& rX,
                                     double u,
                                     Element<DIM,DIM>* pElement=NULL);

    /**
     * Overridden ComputeSourceTermAtNode() method. That is never called.
     *
     * @return computed source term at a node.
     *
     * @param rNode the node at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the node
     */
    virtual double ComputeSourceTermAtNode(const Node<DIM>& rNode, double u);

    /**
     * Overridden ComputeDiffusionTerm() method.
     *
     * @param rX the point in space at which the diffusion term is computed
     * @param pElement the mesh element that x is contained in (optional; defaults to NULL).
     *
     * @return a matrix.
     */
    virtual c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement=NULL);

    /**
     * @return the uptake rate.
     *
     * @param elementIndex the element we wish to return the uptake rate for
     */
    double GetUptakeRateForElement(unsigned elementIndex);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PottsParabolicPde)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a PottsParabolicPde.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const PottsParabolicPde<DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const PottsBasedCellPopulation<DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a PottsParabolicPde.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, PottsParabolicPde<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    PottsBasedCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)PottsParabolicPde<DIM>(*p_cell_population);
}
}
} // namespace ...

#endif /*AVERAGESOURCEPARABOLICPDE_HPP_*/
