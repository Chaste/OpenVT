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

#include "ContactInhibitionModifier_mod.hpp"
#include "FixedDurationCellCycleModel.hpp"

template<unsigned DIM>
ContactInhibitionModifier_mod<DIM>::ContactInhibitionModifier_mod()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
ContactInhibitionModifier_mod<DIM>::~ContactInhibitionModifier_mod()
{
}

template<unsigned DIM>
void ContactInhibitionModifier_mod<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ContactInhibitionModifier_mod<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ContactInhibitionModifier_mod<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // TRACE("UPDATING INHIBITION DATA");
    rCellPopulation.Update(); // Make sure the cell population is updated


    for (typename AbstractCellPopulation<DIM>::Iterator pCell = rCellPopulation.Begin();
         pCell != rCellPopulation.End();
         ++pCell)
    {
        // Determine the deformable radius based on the current radius and the neighbouring cells
        unsigned index = rCellPopulation.GetLocationIndexUsingCell(*pCell);
        std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(index);
        c_vector<double, DIM> cell_position = rCellPopulation.GetLocationOfCellCentre(*pCell);
        double cell_radius = pCell->GetCellData()->GetItem("Radius");

        // Iterate over these neighbours
        double effective_radius_sum = 0.0;
        int num_neighbours_in_contact = 0;
        for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                neighbour_iter != neighbour_indices.end();
                ++neighbour_iter)
        {
            // Get Cell location index
            unsigned neighbour_index = *neighbour_iter;
            if (neighbour_index != index) // Just in case the cell population returns the cell itself as a neighbour
            {
                num_neighbours_in_contact++;
                CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);

                c_vector<double, DIM> neighbouring_cell_position = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell);
                double neighbouring_cell_radius = p_neighbour_cell->GetCellData()->GetItem("Radius");

                // Compute distance between cell centers
                double distance_ij = norm_2(cell_position - neighbouring_cell_position);
                if (distance_ij < (cell_radius + neighbouring_cell_radius)) // Only consider neighbours that are close enough to be in contact
                {
                    effective_radius_sum = effective_radius_sum + 0.5*(cell_radius - neighbouring_cell_radius + distance_ij);
                }
            }
        }
        effective_radius_sum = (1.0/6.0)*(effective_radius_sum + cell_radius*(6 - num_neighbours_in_contact));
        // if (effective_radius_sum > cell_radius)
        // {
        //     PRINT_VARIABLE(num_neighbours_in_contact);
        //     PRINT_2_VARIABLES(effective_radius_sum, cell_radius);
        //     PRINT_VECTOR(cell_position);
            
        //     double effective_radius_sum = 0.0;
        //     for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
        //             neighbour_iter != neighbour_indices.end();
        //             ++neighbour_iter)
        //     {
        //         // Get Cell location index
        //         unsigned neighbour_index = *neighbour_iter;
        //         CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);

        //         c_vector<double, DIM> neighbouring_cell_position = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell);
        //         double neighbouring_cell_radius = p_neighbour_cell->GetCellData()->GetItem("Radius");

        //         // Compute distance between cell centers
        //         double distance_ij = norm_2(cell_position - neighbouring_cell_position);
        //         if (distance_ij < (cell_radius + neighbouring_cell_radius)) // Only consider neighbours that are close enough to be in contact
        //         {
        //             effective_radius_sum = effective_radius_sum + 0.5*(cell_radius - neighbouring_cell_radius + distance_ij);
        //             PRINT_2_VARIABLES(distance_ij, neighbouring_cell_radius);
        //         }
        //     }
        //     effective_radius_sum = (1.0/6.0)*(effective_radius_sum + cell_radius*(6 - neighbour_indices.size()));
            
            
        // }
        pCell->GetCellData()->SetItem("Deformable Radius", effective_radius_sum);
    }
    

    for (typename AbstractCellPopulation<DIM>::Iterator pCell = rCellPopulation.Begin();
         pCell != rCellPopulation.End();
         ++pCell)
    {
        // Get the location index corresponding to this cell
        unsigned index = rCellPopulation.GetLocationIndexUsingCell(*pCell);
        std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(index);

        c_vector<double, DIM> cell_position = rCellPopulation.GetLocationOfCellCentre(*pCell);

        double cell_radius = pCell->GetCellData()->GetItem("Radius");
        // double cell_radius = pCell->GetCellData()->GetItem("Deformable Radius");

        double relative_surface_sum = 0.0;
        double relative_area_sum = 0.0;

        // Iterate over these neighbours
        for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
             neighbour_iter != neighbour_indices.end();
             ++neighbour_iter)
        {
            // Get Cell location index
            unsigned neighbour_index = *neighbour_iter;
            if (neighbour_index != index) // Just in case the cell population returns the cell itself as a neighbour
            {
                CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);

                c_vector<double, DIM> neighbouring_cell_position = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell);

                
                double neighbouring_cell_radius = p_neighbour_cell->GetCellData()->GetItem("Radius");
                // double neighbouring_cell_radius = p_neighbour_cell->GetCellData()->GetItem("Deformable Radius");

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

        pCell->GetCellData()->SetItem("FreeSurfaceFraction", std::max(0.0, free_surface_fraction));
        pCell->GetCellData()->SetItem("FreeAreaFraction",    std::max(0.0, free_area_fraction));




        /********************************************************************************/
        // double free_surface_fraction = 1.0;
        // double free_area_fraction = 0.0;

        // // Iterate over these neighbours
        // for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
        //      neighbour_iter != neighbour_indices.end();
        //      ++neighbour_iter)
        // {
        //     // Get Cell location index
        //     unsigned neighbour_index = *neighbour_iter;
        //     if (neighbour_index != index) // Just in case the cell population returns the cell itself as a neighbour
        //     {
        //         CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);

        //         c_vector<double, DIM> neighbouring_cell_position = rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell);
                
        //         // double neighbouring_cell_radius = p_neighbour_cell->GetCellData()->GetItem("Radius");
        //         double neighbouring_cell_radius = p_neighbour_cell->GetCellData()->GetItem("Deformable Radius");

        //         // Compute distance between cell centers
        //         double distance_ij = norm_2(cell_position - neighbouring_cell_position);

        //         if (distance_ij < (cell_radius + neighbouring_cell_radius)) // Only consider neighbours that are close enough to be in contact
        //         {
        //             double a_i = (1.0/distance_ij) * std::sqrt(4*distance_ij*distance_ij*cell_radius*cell_radius - std::pow(distance_ij*distance_ij - neighbouring_cell_radius*neighbouring_cell_radius + cell_radius*cell_radius, 2));
        //             double theta_i = 2.0*std::asin(a_i /(2.0*cell_radius));
                    
        //             free_surface_fraction = free_surface_fraction - (theta_i / (2.0* M_PI));                   
        //         }
        //     }
        // }
        // free_area_fraction = ((pCell->GetCellData()->GetItem("Deformable Radius"))*(pCell->GetCellData()->GetItem("Deformable Radius"))) / ((pCell->GetCellData()->GetItem("Radius"))*(pCell->GetCellData()->GetItem("Radius")));
        // pCell->GetCellData()->SetItem("FreeSurfaceFraction", std::max(0.0, free_surface_fraction));
        // pCell->GetCellData()->SetItem("FreeAreaFraction",    std::max(0.0, free_area_fraction));

    }

    // TRACE("FINISHED INHIBITION DATA");
    rCellPopulation.Update(); // Make sure the cell population is updated

}

template<unsigned DIM>
void ContactInhibitionModifier_mod<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ContactInhibitionModifier_mod<2>;
template class ContactInhibitionModifier_mod<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ContactInhibitionModifier_mod)
