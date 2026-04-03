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

#include "GeneralisedLinearSpringForceWithMinDistanceItem_mod.hpp"

#include "AbstractCentreBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Debug.hpp"
#include "AbstractCellCycleModel.hpp"
#include <cmath>

#include "FixedDurationCellCycleModelWithContactInhibition.hpp"
#include "FixedGrowthModelWithContactInhibition.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GeneralisedLinearSpringForceWithMinDistanceItem_mod<ELEMENT_DIM,SPACE_DIM>::GeneralisedLinearSpringForceWithMinDistanceItem_mod()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
     mMeinekeSpringStiffness(5.0),        // denoted by mu in Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
     mMeinekeDivisionRestingSpringLength(0.1),
     mMeinekeSpringGrowthDuration(1.0),
     mForceLawType("linear")
{
    if (SPACE_DIM == 1)
    {
        mMeinekeSpringStiffness = 5.0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithMinDistanceItem_mod<ELEMENT_DIM,SPACE_DIM>::VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                                                     unsigned nodeBGlobalIndex,
                                                                                     AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                                                     bool isCloserThanRestLength)
{
    return 1.0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GeneralisedLinearSpringForceWithMinDistanceItem_mod<ELEMENT_DIM,SPACE_DIM>::~GeneralisedLinearSpringForceWithMinDistanceItem_mod()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> GeneralisedLinearSpringForceWithMinDistanceItem_mod<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // TRACE("CALCULATING FORCE BETWEEN NODES");

    // *************************************************************************** //
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<SPACE_DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<SPACE_DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Get the node locations
    const c_vector<double, SPACE_DIM>& r_node_a_location = p_node_a->rGetLocation();
    const c_vector<double, SPACE_DIM>& r_node_b_location = p_node_b->rGetLocation();

    // Get the node radii for a NodeBasedCellPopulation
    double node_a_radius = 0.0;
    double node_b_radius = 0.0;

    // Update actual cell radius
    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    double current_radius_a = p_cell_A->GetCellData()->GetItem("Radius");
    // double current_radius_a = p_cell_A->GetCellData()->GetItem("Deformable Radius");
    p_node_a->SetRadius(current_radius_a);
    double current_radius_b = p_cell_B->GetCellData()->GetItem("Radius");
    // double current_radius_b = p_cell_B->GetCellData()->GetItem("Deformable Radius");
    p_node_b->SetRadius(current_radius_b);
    // TRACE("Collected Deformable Radius from cell data");

    if (bool(dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        // node_a_radius = p_node_a->GetRadius();
        // node_b_radius = p_node_b->GetRadius();
        node_a_radius = current_radius_a;
        node_b_radius = current_radius_b;
    }

    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, SPACE_DIM> unit_difference;
    /*
     * We use the mesh method GetVectorFromAtoB() to compute the direction of the
     * unit vector along the line joining the two nodes, rather than simply subtract
     * their positions, because this method can be overloaded (e.g. to enforce a
     * periodic boundary in Cylindrical2dMesh).
     */
    unit_difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_node_a_location, r_node_b_location);

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    /*
     * If mUseCutOffLength has been set, then there is zero force between
     * two nodes located a distance apart greater than mMechanicsCutOffLength in AbstractTwoBodyInteractionForce.
     */
    if (this->mUseCutOffLength)
    {
        if (distance_between_nodes >= this->GetCutOffLength())
        {
            return zero_vector<double>(SPACE_DIM); // c_vector<double,SPACE_DIM>() is not guaranteed to be fresh memory
        }
    }

    /*
     * Calculate the rest length of the spring connecting the two nodes with a default
     * value of 1.0.
     */
    double rest_length_final = 1.0;

    if (bool(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)))
    {
        rest_length_final = static_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)->GetRestLength(nodeAGlobalIndex, nodeBGlobalIndex);
    }
    else if (bool(dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        assert(node_a_radius > 0 && node_b_radius > 0);
        rest_length_final = node_a_radius+node_b_radius;
        // rest_length_final = p_cell_A->GetCellData()->GetItem("Radius") + p_cell_B->GetCellData()->GetItem("Radius");
    }

    // TRACE("Calculated rest length");

    double rest_length = rest_length_final;    

    // Get the cell cycle phase
    // FixedGrowthModelWithContactInhibition* p_model_A = dynamic_cast<FixedGrowthModelWithContactInhibition*>(p_cell_A->GetCellCycleModel());
    // FixedGrowthModelWithContactInhibition* p_model_B = dynamic_cast<FixedGrowthModelWithContactInhibition*>(p_cell_B->GetCellCycleModel());

    // double phase_timer_A = p_model_A->GetPhaseTimer();


    // double ageA = p_cell_A->GetAge();
    // double ageB = p_cell_B->GetAge();

    // assert(!std::isnan(ageA));
    // assert(!std::isnan(ageB));

    // AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);
    // std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

    // if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
    // {
    //     // Spring rest length increases from a small value to the normal rest length over 1 hour
    //     double lambda = mMeinekeDivisionRestingSpringLength;
    //     rest_length = lambda + (rest_length_final - lambda) * phase_timer_A/mMeinekeSpringGrowthDuration;
    // }

    // if (phase_timer_A + SimulationTime::Instance()->GetTimeStep() >= mMeinekeSpringGrowthDuration && p_static_cast_cell_population->IsMarkedSpring(cell_pair))
    // {
    //     // This spring is about to go out of scope
    //     p_static_cast_cell_population->UnmarkSpring(cell_pair);
    // }

    // TRACE("Done Cell Pair Marking");

    /*
     * If the cells are both newly divided, then the rest length of the spring
     * connecting them grows linearly with time, until 1 hour after division.
     */
    // if (ageA < mMeinekeSpringGrowthDuration && ageB < mMeinekeSpringGrowthDuration)
    // {

    //     // PRINT_3_VARIABLES(ageA,ageB,SimulationTime::Instance()->GetTimeStep());
        
    //     if(ageA == ageB && ageA <= SimulationTime::Instance()->GetTimeStep() )
    //     {
    //         // PRINT_VECTOR(r_node_a_location);
    //         // PRINT_VECTOR(r_node_b_location);
    //         if (norm_2(r_node_a_location-r_node_b_location) < 1.05*mMeinekeDivisionRestingSpringLength)
    //         {
    //             // This spring has just been created
    //             // PRINT_VARIABLE("Marking spring");
    //             p_static_cast_cell_population->MarkSpring(cell_pair);
    //             p_model_A->SetPhaseTimer(0.0);
    //             p_model_B->SetPhaseTimer(0.0);
    //             p_cell_A->GetCellData()->SetItem("cell age", 0.0);
    //             p_cell_B->GetCellData()->SetItem("cell age", 0.0);

    //         }
            
    //     }
        
    // }
    // TRACE("Calculated rest length for newly divided cells");

    /*
     * For apoptosis, progressively reduce the radius of the cell
     */
    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    if (bool(dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        assert(node_a_radius > 0 && node_b_radius > 0);
        a_rest_length = (node_a_radius/(node_a_radius+node_b_radius))*rest_length;
        b_rest_length = (node_b_radius/(node_a_radius+node_b_radius))*rest_length;
    }

    /*
     * If either of the cells has begun apoptosis, then the length of the spring
     * connecting them decreases linearly with time.
     */
    if (p_cell_A->HasApoptosisBegun())
    {
        double time_until_death_a = p_cell_A->GetTimeUntilDeath();
        a_rest_length = a_rest_length * time_until_death_a / p_cell_A->GetApoptosisTime();
    }
    if (p_cell_B->HasApoptosisBegun())
    {
        double time_until_death_b = p_cell_B->GetTimeUntilDeath();
        b_rest_length = b_rest_length * time_until_death_b / p_cell_B->GetApoptosisTime();
    }

    // PRINT_VARIABLE(a_rest_length);
    rest_length = a_rest_length + b_rest_length;
    rest_length_final = rest_length;
    //assert(rest_length <= 1.0+1e-12); ///\todo #1884 Magic number: would "<= 1.0" do?

    // Although in this class the 'spring constant' is a constant parameter, in
    // subclasses it can depend on properties of each of the cells
    double overlap = distance_between_nodes - rest_length;
    bool is_closer_than_rest_length = (overlap <= 0);
    double multiplication_factor = VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex, nodeBGlobalIndex, rCellPopulation, is_closer_than_rest_length);
    double spring_stiffness = mMeinekeSpringStiffness;

    // TRACE("Let's calculate the force between the nodes");

    if (mForceLawType == "linear")
    {
        if (is_closer_than_rest_length)
        {
            return multiplication_factor * spring_stiffness * unit_difference * overlap;        }
        else
        {
            return zero_vector<double>(SPACE_DIM);
        }
        
    }
    else if (mForceLawType == "quadratic")
    {
        if (is_closer_than_rest_length)
        {
            return - multiplication_factor * spring_stiffness * unit_difference * pow((1.0 - distance_between_nodes/rest_length),2.0);
        }
        else
        {
            return zero_vector<double>(SPACE_DIM);
        }
    }
    else if (mForceLawType == "log")
    {
        if (is_closer_than_rest_length)
        {
            assert(overlap > -rest_length_final);
            return multiplication_factor*spring_stiffness * unit_difference * rest_length_final* log(1.0 + overlap/rest_length_final);
        }
        else
        {
            return zero_vector<double>(SPACE_DIM);
        }
    }
    else
    {
        return zero_vector<double>(SPACE_DIM);
    }

    // if (bool(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)))
    // {
    //     return multiplication_factor * spring_stiffness * unit_difference * overlap;
    // }
    // else
    // {
    //     // A reasonably stable simple force law
    //     if (is_closer_than_rest_length) //overlap is negative
    //     {
    //         //log(x+1) is undefined for x<=-1
    //         assert(overlap > -rest_length_final);
    //         c_vector<double, SPACE_DIM> temp = multiplication_factor*spring_stiffness * unit_difference * rest_length_final* log(1.0 + overlap/rest_length_final);
    //         return temp;
    //     }
    //     else
    //     {
    //         // double alpha = 5.0;
    //         // c_vector<double, SPACE_DIM> temp = multiplication_factor*spring_stiffness * unit_difference * overlap * exp(-alpha * overlap/rest_length_final);
    //         // return temp;
    //         return zero_vector<double>(SPACE_DIM);
    //     }
    // }
    // TRACE("FINISHED CALCULATING FORCE BETWEEN NODES");

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithMinDistanceItem_mod<ELEMENT_DIM,SPACE_DIM>::GetMeinekeSpringStiffness()
{
    return mMeinekeSpringStiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithMinDistanceItem_mod<ELEMENT_DIM,SPACE_DIM>::GetMeinekeDivisionRestingSpringLength()
{
    return mMeinekeDivisionRestingSpringLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithMinDistanceItem_mod<ELEMENT_DIM,SPACE_DIM>::GetMeinekeSpringGrowthDuration()
{
    return mMeinekeSpringGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithMinDistanceItem_mod<ELEMENT_DIM,SPACE_DIM>::SetMeinekeSpringStiffness(double springStiffness)
{
    assert(springStiffness > 0.0);
    mMeinekeSpringStiffness = springStiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithMinDistanceItem_mod<ELEMENT_DIM,SPACE_DIM>::SetMeinekeDivisionRestingSpringLength(double divisionRestingSpringLength)
{
    assert(divisionRestingSpringLength <= 1.0);
    assert(divisionRestingSpringLength >= 0.0);

    mMeinekeDivisionRestingSpringLength = divisionRestingSpringLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithMinDistanceItem_mod<ELEMENT_DIM,SPACE_DIM>::SetMeinekeSpringGrowthDuration(double springGrowthDuration)
{
    assert(springGrowthDuration >= 0.0);

    mMeinekeSpringGrowthDuration = springGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithMinDistanceItem_mod<ELEMENT_DIM,SPACE_DIM>::SetForceLawType(std::string force_law)
{
    mForceLawType = force_law;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithMinDistanceItem_mod<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MeinekeSpringStiffness>" << mMeinekeSpringStiffness << "</MeinekeSpringStiffness>\n";
    *rParamsFile << "\t\t\t<MeinekeDivisionRestingSpringLength>" << mMeinekeDivisionRestingSpringLength << "</MeinekeDivisionRestingSpringLength>\n";
    *rParamsFile << "\t\t\t<MeinekeSpringGrowthDuration>" << mMeinekeSpringGrowthDuration << "</MeinekeSpringGrowthDuration>\n";
    *rParamsFile << "\t\t\t<ForceLawType>" << mForceLawType << "</ForceLawType>\n";

    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class GeneralisedLinearSpringForceWithMinDistanceItem_mod<1,1>;
template class GeneralisedLinearSpringForceWithMinDistanceItem_mod<1,2>;
template class GeneralisedLinearSpringForceWithMinDistanceItem_mod<2,2>;
template class GeneralisedLinearSpringForceWithMinDistanceItem_mod<1,3>;
template class GeneralisedLinearSpringForceWithMinDistanceItem_mod<2,3>;
template class GeneralisedLinearSpringForceWithMinDistanceItem_mod<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(GeneralisedLinearSpringForceWithMinDistanceItem_mod)
