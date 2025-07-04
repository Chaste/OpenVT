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

#ifndef ABSTRACTPottsPDEMODIFIER_HPP_
#define ABSTRACTPottsPDEMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractPdeModifier.hpp"
#include "BoundaryConditionsContainer.hpp"

/**
 * An abstract modifier class containing functionality common to EllipticPottsPdeModifier
 * and ParabolicPottsPdeModifier, which both solve a linear elliptic or parabolic PDE
 * coupled to a cell-based simulation on a coarse domain.
 */
template<unsigned DIM>
class AbstractPottsPdeModifier : public AbstractPdeModifier<DIM>
{
    friend class TestEllipticPottsPdeModifier;
    friend class TestParabolicPottsPdeModifier;
    friend class TestOffLatticeSimulationWithPdes;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractPdeModifier<DIM> >(*this);
        archive & mpMeshCuboid;
        archive & mStepSize;
        archive & mSetBcsOnBoxBoundary;
    }

protected:

    /** Map between cells and the elements of the FE mesh containing them. */
    std::map<CellPtr, unsigned> mCellPdeElementMap;

    /**
     * Pointer to a ChasteCuboid storing the outer boundary for the FE mesh.
     * Stored as a member to facilitate archiving.
     */
    boost::shared_ptr<ChasteCuboid<DIM> > mpMeshCuboid;

    /**
     * The step size to be used in the FE mesh.
     * Stored as a member to facilitate archiving.
     */
    double mStepSize;

    /**
     * Whether to set the boundary condition on the edge of the box domain rather than the cell population.
     * Default to true.
     */
    bool mSetBcsOnBoxBoundary;
    
public:

    /**
     * Constructor.
     *
     * @param pPde A shared pointer to a linear PDE object (defaults to NULL)
     * @param pBoundaryCondition A shared pointer to an abstract boundary condition
     *     (defaults to NULL, corresponding to a constant boundary condition with value zero)
     * @param isNeumannBoundaryCondition Whether the boundary condition is Neumann (defaults to true)
     * @param pMeshCuboid A shared pointer to a ChasteCuboid specifying the outer boundary for the FE mesh (defaults to NULL)
     * @param stepSize step size to be used in the FE mesh (defaults to 1.0, i.e. the default cell size)
     * @param solution solution vector (defaults to NULL)
     */
    AbstractPottsPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde=boost::shared_ptr<AbstractLinearPde<DIM,DIM> >(),
                                 boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition=boost::shared_ptr<AbstractBoundaryCondition<DIM> >(),
                                 bool isNeumannBoundaryCondition=true,
                                 boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid=boost::shared_ptr<ChasteCuboid<DIM> >(),
                                 double stepSize=1.0,
                                 Vec solution=nullptr);

    /**
     * Destructor.
     */
    virtual ~AbstractPottsPdeModifier();

    /**
     * @return mStepSize.
     */
    double GetStepSize();

    /**
     * Set mSetBcsOnCoarseBoundary.
     *
     * @param setBcsOnBoxBoundary whether to set the boundary condition on the edge of the box domain rather than the cell population
     */
    void SetBcsOnBoxBoundary(bool setBcsOnBoxBoundary);

    /**
     * @return mSetBcsOnCoarseBoundary.
     */
    bool AreBcsSetOnBoxBoundary();

    /**
     * Helper method to construct the boundary conditions container for the PDE.
     *
     * @param rCellPopulation reference to the cell population
     * @param pBcc the boundary conditions container to fill 
     */
    void ConstructBoundaryConditionsContainerHelper(AbstractCellPopulation<DIM,DIM>& rCellPopulation, 
                                                    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > pBcc);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * Here we just initialize the Cell PDE element map
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to generate the pde mesh for the first time.
     *
     * @param pMeshCuboid the outer boundary for the FE mesh.
     * @param stepSize the step size to be used in the FE mesh.
     */
    void GenerateFeMesh(boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid, double stepSize);

    /**
     * Helper method to generate a pde mesh.
     *
     * @param pMeshCuboid the outer boundary for the FE mesh.
     * @param stepSize the step size to be used in the FE mesh.
     * @param pMesh a pointer to the mesh to be created
     */
    void GenerateAndReturnFeMesh(boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid, double stepSize, TetrahedralMesh<DIM,DIM>* pMesh);

    /**
     * Helper method to copy the PDE solution to CellData
     *
     * Here we need to interpolate from the FE mesh onto the cells.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractPottsPdeModifier)

#endif /*ABSTRACTPottsPDEMODIFIER_HPP_*/
