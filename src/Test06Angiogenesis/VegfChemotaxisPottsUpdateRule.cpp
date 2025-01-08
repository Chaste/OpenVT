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

#include "VegfChemotaxisPottsUpdateRule.hpp"
#include "Debug.hpp"

template<unsigned DIM>
VegfChemotaxisPottsUpdateRule<DIM>::VegfChemotaxisPottsUpdateRule()
    : AbstractPottsUpdateRule<DIM>(),
      mChemotaxisParameter(10.0) // Educated guess
{
}

template<unsigned DIM>
VegfChemotaxisPottsUpdateRule<DIM>::~VegfChemotaxisPottsUpdateRule()
{
}

template<unsigned DIM>
double VegfChemotaxisPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                        unsigned targetNodeIndex,
                                                                        PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    std::set<unsigned> containing_elements = rCellPopulation.GetNode(currentNodeIndex)->rGetContainingElementIndices();
    std::set<unsigned> new_location_containing_elements = rCellPopulation.GetNode(targetNodeIndex)->rGetContainingElementIndices();

    // Every node must each be in at most one element
    assert(containing_elements.size() < 2);
    assert(new_location_containing_elements.size() < 2);


    bool current_node_contained = !containing_elements.empty();
    bool target_node_contained = !new_location_containing_elements.empty();

    if (!current_node_contained && !target_node_contained)
    {
        EXCEPTION("At least one of the current node or target node must be in an element.");
    }

    double delta_H = 0.0;

    if (current_node_contained && target_node_contained)
    {
        // both cells so return 0 
    }
    else 
    {
        assert((!current_node_contained && target_node_contained) || (current_node_contained &&!target_node_contained));
        
        double current_vegf = rCellPopulation.GetPdeSolutionAtNode(currentNodeIndex);

        double target_vegf = rCellPopulation.GetPdeSolutionAtNode(targetNodeIndex);
//PRINT_2_VARIABLES(current_vegf,target_vegf);

        delta_H = mChemotaxisParameter * (current_vegf - target_vegf);

    }
    return delta_H;
}

template<unsigned DIM>
double VegfChemotaxisPottsUpdateRule<DIM>::GetChemotaxisParameter()
{
    return mChemotaxisParameter;
}

template<unsigned DIM>
void VegfChemotaxisPottsUpdateRule<DIM>::SetChemotaxisParameter(double chemotaxisParameter)
{
    mChemotaxisParameter = chemotaxisParameter;
}


template<unsigned DIM>
void VegfChemotaxisPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ChemotaxisParameter>" << mChemotaxisParameter << "</ChemotaxisParameter>\n";

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

// Explicit instantiation
template class VegfChemotaxisPottsUpdateRule<1>;
template class VegfChemotaxisPottsUpdateRule<2>;
template class VegfChemotaxisPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VegfChemotaxisPottsUpdateRule)
