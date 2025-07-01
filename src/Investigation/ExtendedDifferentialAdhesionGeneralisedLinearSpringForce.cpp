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

#include "ExtendedDifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "CellLabel.hpp"
#include "MediumCellProperty.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::ExtendedDifferentialAdhesionGeneralisedLinearSpringForce()
   : GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>(),
     mHomotypicLabelledSpringConstantMultiplier(1.0),
     mHeterotypicSpringConstantMultiplier(1.0),
     mLabelledMediumSpringConstantMultiplier(1.0),
     mUnLabelledMediumSpringConstantMultiplier(1.0),
     mMediumMediumSpringConstantMultiplier(1.0)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::VariableSpringConstantMultiplicationFactor(
    unsigned nodeAGlobalIndex,
    unsigned nodeBGlobalIndex,
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,
    bool isCloserThanRestLength)
{

    if (isCloserThanRestLength)
    {
        return 1.0;
    }
    else
    {
        // Determine which (if any) of the cells corresponding to these nodes are labelled
        CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
        bool cell_A_is_labelled = p_cell_A->template HasCellProperty<CellLabel>();

        CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);
        bool cell_B_is_labelled = p_cell_B->template HasCellProperty<CellLabel>();

        // Determine which (if any) of the cells corresponding to these nodes are medium
        bool cell_A_is_medium = p_cell_A->template HasCellProperty<MediumCellProperty>();
        bool cell_B_is_medium = p_cell_B->template HasCellProperty<MediumCellProperty>();

        if (cell_A_is_medium && cell_B_is_medium)
        {
            return mMediumMediumSpringConstantMultiplier;
        }
        else if (cell_A_is_medium && !cell_B_is_medium)
        {
            if (cell_B_is_labelled)
            {
                return mLabelledMediumSpringConstantMultiplier;
            }
            else
            {
                return mUnLabelledMediumSpringConstantMultiplier;
            }
        }
        else if (!cell_A_is_medium && cell_B_is_medium)
        {
            if (cell_A_is_labelled)
            {
                return mLabelledMediumSpringConstantMultiplier;
            }
            else
            {
                return mUnLabelledMediumSpringConstantMultiplier;
            }
        }
        else //(!cell_A_is_medium && !cell_B_is_medium)
        {   // For heterotypic interactions between  cells, scale the spring constant by mHeterotypicSpringConstantMultiplier
            if (cell_A_is_labelled != cell_B_is_labelled)
            {
                return mHeterotypicSpringConstantMultiplier;
            }
            else
            {
                // For homotypic interactions between 2 labelled cells, scale the spring constant by mHomotypicLabelledSpringConstantMultiplier
                if (cell_A_is_labelled)
                {
                    return mHomotypicLabelledSpringConstantMultiplier;
                }
                else
                {
                    // For homotypic interactions between unlabelled cells, leave the spring constant unchanged from its normal value
                    return 1.0;
                }
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::GetHomotypicLabelledSpringConstantMultiplier()
{
    return mHomotypicLabelledSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::SetHomotypicLabelledSpringConstantMultiplier(double labelledSpringConstantMultiplier)
{
    assert(labelledSpringConstantMultiplier > 0.0);
    mHomotypicLabelledSpringConstantMultiplier = labelledSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::GetHeterotypicSpringConstantMultiplier()
{
    return mHeterotypicSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::SetHeterotypicSpringConstantMultiplier(double heterotypicSpringConstantMultiplier)
{
    assert(springConstantMultiplier > 0.0);
    mHeterotypicSpringConstantMultiplier = heterotypicSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::GetLabelledMediumSpringConstantMultiplier()
{
    return mLabelledMediumSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::SetLabelledMediumSpringConstantMultiplier(double labelledMediumSpringConstantMultiplier)
{
    assert(labelledMediumSpringConstantMultiplier > 0.0);
    mLabelledMediumSpringConstantMultiplier = labelledMediumSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::GetUnLabelledMediumSpringConstantMultiplier()
{
    return mUnLabelledMediumSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::SetUnLabelledMediumSpringConstantMultiplier(double unLabelledMediumSpringConstantMultiplier)
{
    assert(unLabelledMediumSpringConstantMultiplier > 0.0);
    mUnLabelledMediumSpringConstantMultiplier = unLabelledMediumSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::GetMediumMediumSpringConstantMultiplier()
{
    return mMediumMediumSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::SetMediumMediumSpringConstantMultiplier(double mediumMediumSpringConstantMultiplier)
{
    assert(mediumMediumSpringConstantMultiplier > 0.0);
    mMediumMediumSpringConstantMultiplier = mediumMediumSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<HomotypicLabelledSpringConstantMultiplier>" << mHomotypicLabelledSpringConstantMultiplier << "</HomotypicLabelledSpringConstantMultiplier>\n";
    *rParamsFile << "\t\t\t<HeterotypicSpringConstantMultiplier>" << mHeterotypicSpringConstantMultiplier << "</HeterotypicSpringConstantMultiplier>\n";

    // Call direct parent class
    GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<1,1>;
template class ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<1,2>;
template class ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<2,2>;
template class ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<1,3>;
template class ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<2,3>;
template class ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ExtendedDifferentialAdhesionGeneralisedLinearSpringForce)
