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

template<unsigned DIM>
PottsParabolicPde<DIM>::PottsParabolicPde(PottsBasedCellPopulation<DIM>& rCellPopulation,
                                                            double constantCellSourceCoefficient, 
                                                            double linearCellSourceCoefficient, 
                                                            double constantSourceCoefficient, 
                                                            double linearSourceCoefficient, 
                                                            double diffusionCoefficient,
                                                            double duDtCoefficient,
                                                            bool scaleByCellVolume)
    : mrCellPopulation(rCellPopulation),
      mConstantCellSourceCoefficient(constantCellSourceCoefficient),
      mLinearCellSourceCoefficient(linearCellSourceCoefficient),
      mConstantSourceCoefficient(constantSourceCoefficient),
      mLinearSourceCoefficient(linearSourceCoefficient),
      mDiffusionCoefficient(diffusionCoefficient),
      mDuDtCoefficient(duDtCoefficient),
      mScaleByCellVolume(scaleByCellVolume)
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
bool PottsParabolicPde<DIM>::GetScaleByCellVolume() const
{
    return mScaleByCellVolume;
}

template<unsigned DIM>
void PottsParabolicPde<DIM>::SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap) // must be called before solve
{
    NEVER_REACHED
    // Allocate memory
    mCellDensityOnCoarseElements.resize(rCoarseMesh.GetNumElements());
    for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
    {
        mCellDensityOnCoarseElements[elem_index] = 0.0;
    }

    {
        // Loop over cells, find which coarse element it is in, and add 1 to mSourceTermOnCoarseElements[elem_index]
        for (typename PottsBasedCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
            cell_iter != mrCellPopulation.End();
            ++cell_iter)
        {
            unsigned elem_index = UNSIGNED_UNSET;
            const ChastePoint<DIM>& r_position_of_cell = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

            if (pCellPdeElementMap != nullptr)
            {
                elem_index = (*pCellPdeElementMap)[*cell_iter];
            }
            else
            {
                elem_index = rCoarseMesh.GetContainingElementIndex(r_position_of_cell);
            }

            bool cell_is_apoptotic = cell_iter->template HasCellProperty<ApoptoticCellProperty>();

            if (!cell_is_apoptotic)
            {
                double cell_weight = 1.0;

                if (mScaleByCellVolume)
                {   
                    // If scaling by cell volume then use volume here instead of cell count 
                    cell_weight = mrCellPopulation.GetVolumeOfCell(*cell_iter);

                    if (cell_weight <1e-6)
                    {
                        EXCEPTION("The volume of one of the cells is " << cell_weight << 
                                " and you are scaling by cell volume. Either turn scaling off or use"  
                                " a cell model with non zero areas (i.e. a Bounded Voronoi Tesselation model).");
                    }
                }

                mCellDensityOnCoarseElements[elem_index] += cell_weight;
            }
        }
    }

    // Then divide each entry of mSourceTermOnCoarseElements by the element's area
    c_matrix<double, DIM, DIM> jacobian;
    double det;
    for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
    {
        rCoarseMesh.GetElement(elem_index)->CalculateJacobian(jacobian, det);
        mCellDensityOnCoarseElements[elem_index] /= rCoarseMesh.GetElement(elem_index)->GetVolume(det);
    }
}

template<unsigned DIM>
double PottsParabolicPde<DIM>::ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& )
{
    return mDuDtCoefficient;
}

template<unsigned DIM>
double PottsParabolicPde<DIM>::ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement)
{
    
    assert(!mCellDensityOnCoarseElements.empty());

    // The source term is (a*density + c)*u + b*density + d
    bool is_in_cell = 1;

    double constant_source_term = mConstantCellSourceCoefficient * is_in_cell + mConstantSourceCoefficient;
    double linear_source_term_coeficient = mLinearCellSourceCoefficient * is_in_cell + mLinearSourceCoefficient;
    
    double source_term = constant_source_term + linear_source_term_coeficient * u;

    return source_term;
    //return 25*mCellDensityOnCoarseElements[pElement->GetIndex()] - u;
}

// LCOV_EXCL_START
template<unsigned DIM>
double PottsParabolicPde<DIM>::ComputeSourceTermAtNode(const Node<DIM>& rNode, double u)
{
    NEVER_REACHED;
    return 0.0;
}
// LCOV_EXCL_STOP

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
