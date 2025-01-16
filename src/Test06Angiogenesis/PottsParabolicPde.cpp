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

#include "PottsParabolicPde.hpp"
#include "ApoptoticCellProperty.hpp"
#include "Debug.hpp"

template<unsigned DIM>
PottsParabolicPde<DIM>::PottsParabolicPde(PottsBasedCellPopulation<DIM>& rCellPopulation,
                                                            double constantCellSourceCoefficient, 
                                                            double linearCellSourceCoefficient, 
                                                            double constantSourceCoefficient, 
                                                            double linearSourceCoefficient, 
                                                            double diffusionCoefficient,
                                                            double duDtCoefficient)
    : mrCellPopulation(rCellPopulation),
      mConstantCellSourceCoefficient(constantCellSourceCoefficient),
      mLinearCellSourceCoefficient(linearCellSourceCoefficient),
      mConstantSourceCoefficient(constantSourceCoefficient),
      mLinearSourceCoefficient(linearSourceCoefficient),
      mDiffusionCoefficient(diffusionCoefficient),
      mDuDtCoefficient(duDtCoefficient)
{
}

template<unsigned DIM>
const PottsBasedCellPopulation<DIM>& PottsParabolicPde<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
double PottsParabolicPde<DIM>::GetConstantCellCoefficient() const
{
    return mConstantCellSourceCoefficient;
}

template<unsigned DIM>
double PottsParabolicPde<DIM>::GetLinearCellCoefficient() const
{
    return mLinearCellSourceCoefficient;
}

template<unsigned DIM>
double PottsParabolicPde<DIM>::GetConstantCoefficient() const
{
    return mConstantSourceCoefficient;
}

template<unsigned DIM>
double PottsParabolicPde<DIM>::GetLinearCoefficient() const
{
    return mLinearSourceCoefficient;
}

template<unsigned DIM>
double PottsParabolicPde<DIM>::GetDiffusionCoefficient() const
{
    return mDiffusionCoefficient;
}

template<unsigned DIM>
double PottsParabolicPde<DIM>::GetDuDtCoefficient() const
{
    return mDuDtCoefficient;
}

template<unsigned DIM>
void PottsParabolicPde<DIM>::SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap) // must be called before solve
{
    NEVER_REACHED;
}

template<unsigned DIM>
double PottsParabolicPde<DIM>::ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& )
{
    return mDuDtCoefficient;
}

template<unsigned DIM>
double PottsParabolicPde<DIM>::ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement)
{
    NEVER_REACHED;
}

template<unsigned DIM>
double PottsParabolicPde<DIM>::ComputeSourceTermAtNode(const Node<DIM>& rNode, double u)
{
    //mrCellPopulation->GetNumNodes();

    // Need to get index as nodes in FE Mesh and PottsMesh are different (but in same location).
    unsigned node_index = rNode.GetIndex();
    
    std::set<unsigned> containing_elements = mrCellPopulation.GetNode(node_index)->rGetContainingElementIndices();
   
    unsigned node_contained = !containing_elements.empty();

//PRINT_2_VARIABLES(node_index,node_contained);

    double constant_source_term = mConstantCellSourceCoefficient * (double) node_contained + mConstantSourceCoefficient;
    double linear_source_term_coeficient = mLinearCellSourceCoefficient * (double) node_contained + mLinearSourceCoefficient;
    
    double source_term = constant_source_term + linear_source_term_coeficient * u;
//PRINT_3_VARIABLES(constant_source_term,linear_source_term_coeficient, source_term);
    return source_term;
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> PottsParabolicPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return mDiffusionCoefficient*identity_matrix<double>(DIM);
}

template<unsigned DIM>
double PottsParabolicPde<DIM>::GetUptakeRateForElement(unsigned elementIndex)
{
    NEVER_REACHED;
    return this->mCellDensityOnCoarseElements[elementIndex];
}

// Explicit instantiation
template class PottsParabolicPde<1>;
template class PottsParabolicPde<2>;
template class PottsParabolicPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PottsParabolicPde)
