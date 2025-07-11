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

#include "RadialDifferentiationModifier.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

template<unsigned DIM>
RadialDifferentiationModifier<DIM>::RadialDifferentiationModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mRadius(DOUBLE_UNSET)
{
}

template<unsigned DIM>
RadialDifferentiationModifier<DIM>::~RadialDifferentiationModifier()
{
}

template<unsigned DIM>
void RadialDifferentiationModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    DifferentiateCells(rCellPopulation);
}

template<unsigned DIM>
void RadialDifferentiationModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    DifferentiateCells(rCellPopulation);
}

template<unsigned DIM>
void RadialDifferentiationModifier<DIM>::DifferentiateCells(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    assert(mRadius!=DOUBLE_UNSET);

    // Iterate over cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get distance from centre of cell population
        double r = norm_2(rCellPopulation.GetLocationOfCellCentre(*cell_iter));

        if (r > mRadius)
        {
            boost::shared_ptr<AbstractCellProperty> p_diff_type =
                    cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->template Get<DifferentiatedCellProliferativeType>();
            cell_iter->SetCellProliferativeType(p_diff_type);
        }
    }
}

template<unsigned DIM>
double RadialDifferentiationModifier<DIM>::GetRadius()
{
    return mRadius;
}

template<unsigned DIM>
void RadialDifferentiationModifier<DIM>::SetRadius(double radius)
{
    mRadius = radius;
}

template<unsigned DIM>
void RadialDifferentiationModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<mRadius>" << mRadius << "</mRadius>\n";

    // Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class RadialDifferentiationModifier<1>;
template class RadialDifferentiationModifier<2>;
template class RadialDifferentiationModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RadialDifferentiationModifier)

