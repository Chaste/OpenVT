/*

Copyright (c) 2005-2024, University of Oxford.
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

#include "ExtendedAdhesionPottsUpdateRule.hpp"
#include "Debug.hpp"

template<unsigned DIM>
ExtendedAdhesionPottsUpdateRule<DIM>::ExtendedAdhesionPottsUpdateRule()
    : AdhesionPottsUpdateRule<DIM>(),
      mCellBoundaryNodeAdhesionEnergyParameter(100000.0) // Educated guess so cells avoid boundary nodes
{
}

template<unsigned DIM>
ExtendedAdhesionPottsUpdateRule<DIM>::~ExtendedAdhesionPottsUpdateRule()
{
}

template<unsigned DIM>
double ExtendedAdhesionPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                unsigned targetNodeIndex,
                                                                PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    std::set<unsigned> containing_elements = rCellPopulation.GetNode(currentNodeIndex)->rGetContainingElementIndices();
    std::set<unsigned> new_location_containing_elements = rCellPopulation.GetNode(targetNodeIndex)->rGetContainingElementIndices();

    bool current_node_contained = !containing_elements.empty();
    bool target_node_contained = !new_location_containing_elements.empty();

    //bool current_node_is_boundary_node = rCellPopulation.GetNode(currentNodeIndex)->IsBoundaryNode();
    bool target_node_is_boundary_node = rCellPopulation.GetNode(targetNodeIndex)->IsBoundaryNode();

    // Every node must each be in at most one element
    assert(new_location_containing_elements.size() < 2);

    if (!current_node_contained && !target_node_contained)
    {
        EXCEPTION("At least one of the current node or target node must be in an element.");
    }

    if (current_node_contained && target_node_contained)
    {
        if (*(new_location_containing_elements.begin()) == *(containing_elements.begin()))
        {
            EXCEPTION("The current node and target node must not be in the same element.");
        }
    }

    // Iterate over nodes neighbouring the target node to work out the contact energy contribution
    double delta_H = 0.0;
    std::set<unsigned> target_neighbouring_node_indices = rCellPopulation.rGetMesh().GetVonNeumannNeighbouringNodeIndices(targetNodeIndex);
    for (std::set<unsigned>::iterator iter = target_neighbouring_node_indices.begin();
         iter != target_neighbouring_node_indices.end();
         ++iter)
    {
        std::set<unsigned> neighbouring_node_containing_elements = rCellPopulation.rGetMesh().GetNode(*iter)->rGetContainingElementIndices();

        // Every node must each be in at most one element
        assert(neighbouring_node_containing_elements.size() < 2);

        bool neighbouring_node_contained = !neighbouring_node_containing_elements.empty();

        /**
         * Before the move, we have a negative contribution (H_0) to the Hamiltonian if:
         * the target node and neighbouring node are NOT contained in the same Potts element;
         * the neighbouring node is contained in a Potts element, but the target node is not; or
         * the target node is contained in a Potts element, but the neighbouring node is not.
         */
        if (neighbouring_node_contained && target_node_contained)
        {
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            unsigned target_element = (*new_location_containing_elements.begin());
            if (target_element != neighbour_element)
            {
                // The nodes are currently contained in different elements
                delta_H -= this->GetCellCellAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(target_element), rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
            }
        }
        else if (neighbouring_node_contained && !target_node_contained)
        {
            // The neighbouring node is contained in a Potts element, but the target node is not
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            delta_H -= this->GetCellBoundaryAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
        }
        else if (!neighbouring_node_contained && target_node_contained)
        {
            // The target node is contained in a Potts element, but the neighbouring node is not
            unsigned target_element = (*new_location_containing_elements.begin());
            delta_H -= this->GetCellBoundaryAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(target_element));
        }

        // First check if the target node is a boundary node as if it is thid is the main iteraction.
        if (target_node_is_boundary_node)
        {
            delta_H += 10000000;//mCellBoundaryNodeAdhesionEnergyParameter;
        }
        else
        {
            /**
             * After the move, we have a positive contribution (H_1) to the Hamiltonian if:
             * the current node and neighbouring node are contained in different Potts elements;
             * the neighbouring node is contained in a Potts element, but the current node is not; or
             * the current node is contained in a Potts element, but the neighbouring node is not.
             */
            if (neighbouring_node_contained && current_node_contained)
            {
                unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
                unsigned current_element = (*containing_elements.begin());
                if (current_element != neighbour_element)
                {
                    // The nodes are currently contained in different elements
                    delta_H += this->GetCellCellAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(current_element),rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
                }
            }
            else if (neighbouring_node_contained && !current_node_contained)
            {
                // The neighbouring node is contained in a Potts element, but the current node is not
                unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
                delta_H += this->GetCellBoundaryAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
            }
            else if (!neighbouring_node_contained && current_node_contained)
            {
                // The current node is contained in a Potts element, but the neighbouring node is not
                unsigned current_element = (*containing_elements.begin());
                delta_H += this->GetCellBoundaryAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(current_element));
            }
        }
    }

    return delta_H;
}


template<unsigned DIM>
double ExtendedAdhesionPottsUpdateRule<DIM>::GetCellBoundaryNodeAdhesionEnergyParameter()
{
    return mCellBoundaryNodeAdhesionEnergyParameter;
}

template<unsigned DIM>
void ExtendedAdhesionPottsUpdateRule<DIM>::SetCellBoundaryNodeAdhesionEnergyParameter(double cellBoundaryNodeAdhesionEnergyParameter)
{
    mCellBoundaryNodeAdhesionEnergyParameter = cellBoundaryNodeAdhesionEnergyParameter;
}

template<unsigned DIM>
void ExtendedAdhesionPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellBoundaryNodeAdhesionEnergyParameter>" << mCellBoundaryNodeAdhesionEnergyParameter << "</CellBoundaryNodeAdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AdhesionPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

// Explicit instantiation
template class ExtendedAdhesionPottsUpdateRule<1>;
template class ExtendedAdhesionPottsUpdateRule<2>;
template class ExtendedAdhesionPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExtendedAdhesionPottsUpdateRule)
