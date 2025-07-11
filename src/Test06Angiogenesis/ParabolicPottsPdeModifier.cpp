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

#include "ParabolicPottsPdeModifier.hpp"
#include "CellBasedParabolicPdeSolver.hpp"
#include "VtkMeshWriter.hpp"
#include "MutableMesh.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "Debug.hpp"

template<unsigned DIM>
ParabolicPottsPdeModifier<DIM>::ParabolicPottsPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde,
                                                                  boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition,
                                                                  bool isNeumannBoundaryCondition,
                                                                  boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid,
                                                                  double stepSize,
                                                                  Vec solution)
    : AbstractPottsPdeModifier<DIM>(pPde,
                                        pBoundaryCondition,
                                        isNeumannBoundaryCondition,
                                        pMeshCuboid,
                                        stepSize,
                                        solution)
{
}

template<unsigned DIM>
ParabolicPottsPdeModifier<DIM>::~ParabolicPottsPdeModifier()
{
}

template<unsigned DIM>
void ParabolicPottsPdeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Set up boundary conditions
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer(rCellPopulation);

    // Use CellBasedParabolicPdeSolver as Potts PDE
    CellBasedParabolicPdeSolver<DIM> solver(this->mpFeMesh,
                                            boost::static_pointer_cast<AbstractLinearParabolicPde<DIM,DIM> >(this->mpPde).get(),
                                            p_bcc.get());

    ///\todo Investigate more than one PDE time step per spatial step
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();
    double dt = p_simulation_time->GetTimeStep();
    solver.SetTimes(current_time,current_time + dt);
    solver.SetTimeStep(dt);

    // Use previous solution as the initial condition
    Vec previous_solution = this->mSolution;
    solver.SetInitialCondition(previous_solution);

    // Note that the linear solver creates a vector, so we have to keep a handle on the old one
    // in order to destroy it
    this->mSolution = solver.Solve();
    PetscTools::Destroy(previous_solution);

    // Now copy solution to cells
    this->UpdateCellData(rCellPopulation);

    // Using a PottsBasedCellPopulation so copy the PDE solution to the lattice sites.
    assert(dynamic_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation));

    // confirm meshes are the same 
    if (this->mpFeMesh->GetNumNodes() == rCellPopulation.GetNumNodes())
    {
        // Store the PDE solution in an accessible form
        ReplicatableVector solution_repl(this->mSolution);

        unsigned num_nodes = rCellPopulation.GetNumNodes();

        std::vector<double> solution_on_lattice_sites(num_nodes, 0);

        for (unsigned node_index = 0u; node_index < num_nodes; node_index++)
        {
            solution_on_lattice_sites[node_index] = solution_repl[node_index];
        }

        // Store the solution in the PottsBasedCellPopulation
        dynamic_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation)->SetPdeSolution(solution_on_lattice_sites);
    }
    else
    {
        PRINT_2_VARIABLES(this->mpFeMesh->GetNumNodes(), rCellPopulation.GetNumNodes());
        NEVER_REACHED;
    }   
}

template<unsigned DIM>
void ParabolicPottsPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractPottsPdeModifier<DIM>::SetupSolve(rCellPopulation,outputDirectory);

    // Copy the cell data to mSolution (this is the initial condition)
    SetupInitialSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > ParabolicPottsPdeModifier<DIM>::ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    this->ConstructBoundaryConditionsContainerHelper(rCellPopulation,p_bcc);
   
    return p_bcc;
}

template<unsigned DIM>
void ParabolicPottsPdeModifier<DIM>::SetupInitialSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Specify homogeneous initial conditions based upon the values stored in CellData.
    // Note need all the CellDataValues to be the same.

    double initial_condition = rCellPopulation.Begin()->GetCellData()->GetItem(this->mDependentVariableName);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        double initial_condition_at_cell = cell_iter->GetCellData()->GetItem(this->mDependentVariableName);
        UNUSED_OPT(initial_condition_at_cell);
        assert(fabs(initial_condition_at_cell - initial_condition)<1e-12);
    }

    // Initialise mSolution
    this->mSolution = PetscTools::CreateAndSetVec(this->mpFeMesh->GetNumNodes(), initial_condition);
}

template<unsigned DIM>
void ParabolicPottsPdeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractPottsPdeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ParabolicPottsPdeModifier<1>;
template class ParabolicPottsPdeModifier<2>;
template class ParabolicPottsPdeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ParabolicPottsPdeModifier)
