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

#include <cmath>

#include "Debug.hpp"
#include "OffLatticeSimulation.hpp"

#include "ContactInhibitionModifier.hpp"
#include "FixedDurationCellCycleModel.hpp"

template<unsigned DIM>
ContactInhibitionModifier<DIM>::ContactInhibitionModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
ContactInhibitionModifier<DIM>::~ContactInhibitionModifier()
{
}

template<unsigned DIM>
void ContactInhibitionModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ContactInhibitionModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ContactInhibitionModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    rCellPopulation.Update(); // Make sure the cell population is updated
    
    for (typename AbstractCellPopulation<DIM>::Iterator pCell = rCellPopulation.Begin();
         pCell != rCellPopulation.End();
         ++pCell)
    {
        // Get the location index corresponding to this cell
        unsigned index = rCellPopulation.GetLocationIndexUsingCell(*pCell);
        std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(index);

        c_vector<double, DIM> cell_position = rCellPopulation.GetLocationOfCellCentre(*pCell);

        // double ageA = pCell->GetAge();
        
        // double cellVolume = rCellPopulation.GetVolumeOfCell(*pCell);
        double cell_radius = pCell->GetCellData()->GetItem("Radius");
        // double cell_radius = rCellPopulation.GetNode(index)->GetRadius();

        double relative_surface_sum = 0.0;
        double relative_area_sum = 0.0;

        // Iterate over these neighbours
        for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
             neighbour_iter != neighbour_indices.end();
             ++neighbour_iter)
        {
            // Get Cell location index
            unsigned neighbour_index = *neighbour_iter;
            if(neighbour_index != index) // Don't consider the cell itself as a neighbour
            {
                CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);

                c_vector<double, DIM> neighbouring_cell_position = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell);
                // double neighbour_cellVolume = rCellPopulation.GetVolumeOfCell(p_neighbour_cell);
                // double neighbouring_cell_radius = std::cbrt((3.0 * neighbour_cellVolume) / (4.0 * M_PI));
                double neighbouring_cell_radius = p_neighbour_cell->GetCellData()->GetItem("Radius");

                // Compute distance between cell centers
                double distance_ij = norm_2(cell_position - neighbouring_cell_position);

                if (distance_ij < (cell_radius + neighbouring_cell_radius)) // Only consider neighbours that are close enough to be in contact
                {
                    double relative_radius = (distance_ij*distance_ij - neighbouring_cell_radius*neighbouring_cell_radius + cell_radius*cell_radius) / (2.0*distance_ij*cell_radius);

                    relative_area_sum += (std::acos(relative_radius) - relative_radius*std::sqrt(1.0 - relative_radius*relative_radius));
                    relative_surface_sum += std::sqrt(1.0 - relative_radius*relative_radius);
                }
            }
        }

        double free_surface_fraction = 1.0 - (relative_surface_sum / M_PI);
        if (free_surface_fraction < 0.0)
        {
            free_surface_fraction = 0.0;
        }
        double free_area_fraction    = 1.0 - (relative_area_sum / M_PI);
        if (free_area_fraction < 0.0)
        {
            free_area_fraction = 0.0;
        }

        // PRINT_2_VARIABLES(free_surface_fraction, free_area_fraction)

        // pCell->GetCellData()->SetItem("FreeSurfaceFraction", std::max(0.0, free_surface_fraction));
        // pCell->GetCellData()->SetItem("FreeAreaFraction",    std::max(0.0, free_area_fraction));
        pCell->GetCellData()->SetItem("FreeSurfaceFraction", free_surface_fraction);
        pCell->GetCellData()->SetItem("FreeAreaFraction",    free_area_fraction);


        // auto pCellCycleModel = static_cast<FixedDurationCellCycleModel*>(pCell->GetCellCycleModel());
        // double phaseG1Duration = pCellCycleModel->GetG1Duration() * 60.0;

        // // Compute target relative volume %
        // double targetRelativeVolume {0.0};

        // // Get current cell age in minutes
        // double cellAge = pCell->GetCellData()->GetItem("cell age") * 60.0;
        // if (cellAge < 0){
        //     cellAge = 0;
        // }

        // if (cellAge < phaseG1Duration)
        // {
        //     targetRelativeVolume = 100.0;

        // } else {
        //     double age = cellAge - phaseG1Duration;
        //     targetRelativeVolume = 100.0 + 0.2850 * age - 0.0002 * age * age;
        // }

        // // Compute target radius
        // double initialRadius = 0.5; // todo: fix magic number from test setup
        // double initialVolume = (4.0 * M_PI * initialRadius * initialRadius * initialRadius) / 3.0;
        // double targetVolume = (targetRelativeVolume * initialVolume) / 100.0;
        // double targetRadius = std::cbrt((3.0 * targetVolume) / (4.0 * M_PI));

        // // Set target radius
        // pCell->GetCellData()->SetItem("Radius", targetRadius);
        // pCell->GetCellData()->SetItem("TargetVolume", targetVolume);

        // // double cellVolume = rCellPopulation.GetVolumeOfCell(*pCell);
        // cell_radius = std::cbrt((3.0 * cellVolume) / (4.0 * M_PI));
        // pCell->GetCellData()->SetItem("Current Radius", cell_radius);
    }

    rCellPopulation.Update(); // Make sure the cell population is updated

}

template<unsigned DIM>
void ContactInhibitionModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ContactInhibitionModifier<2>;
template class ContactInhibitionModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ContactInhibitionModifier)
