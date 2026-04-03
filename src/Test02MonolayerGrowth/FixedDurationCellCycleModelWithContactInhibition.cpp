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

#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"

#include "FixedDurationCellCycleModelWithContactInhibition.hpp"
#include "Debug.hpp"

#include <iostream>

#include "RandomNumberGenerator.hpp"

template<class Archive>
void FixedDurationCellCycleModelWithContactInhibition::serialize(Archive & archive, const unsigned int version)
{
    // Archive cell-cycle model using serialization code from AbstractSimplePhaseBasedCellCycleModel
    archive & boost::serialization::base_object<AbstractSimplePhaseBasedCellCycleModel>(*this);
}

FixedDurationCellCycleModelWithContactInhibition::FixedDurationCellCycleModelWithContactInhibition()
    : AbstractSimplePhaseBasedCellCycleModel()
{
    SetPhaseDurations();
    mPhaseTimer = 0.0;
    mLastCellAge = 0.0;
    // mBeta = 1.0;  // Default threshold for area fraction below which growth is inhibited
    // mGamma = 1.0; // Default threshold for surface area fraction below which growth is inhibited
}

void FixedDurationCellCycleModelWithContactInhibition::ResetForDivision()
{
    AbstractCellCycleModel::ResetForDivision();
    mCurrentCellCyclePhase = G_ONE_PHASE;
    mPhaseTimer = 0.0;
    mLastCellAge = 0.0;
}

void FixedDurationCellCycleModelWithContactInhibition::SetPhaseDurations()
{
    SetStemCellG1Duration(2.0);  
    SetTransitCellG1Duration(2.0);
    SetSDuration(3.0);
    SetG2Duration(3.0);
    SetMDuration(2.0);
    SetMinimumGapDuration(3.0);
}

void FixedDurationCellCycleModelWithContactInhibition::SetG1Duration()
{
    assert(mpCell != NULL);  // Make sure cell exists

    mG1Duration = 7.0;
}

AbstractCellCycleModel* FixedDurationCellCycleModelWithContactInhibition::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    FixedDurationCellCycleModelWithContactInhibition* pCellCycleModel = new FixedDurationCellCycleModelWithContactInhibition();
    return pCellCycleModel;
}

void FixedDurationCellCycleModelWithContactInhibition::SetPhaseTimer(const double value) {
    mPhaseTimer = value;
}

double FixedDurationCellCycleModelWithContactInhibition::GetPhaseTimer()
{
    return mPhaseTimer;
}

bool FixedDurationCellCycleModelWithContactInhibition::ReadyToDivide() {
    assert(mpCell != nullptr);
    
    if (!mReadyToDivide)
    {
        UpdateCellCyclePhase();
        if ((mCurrentCellCyclePhase != G_ZERO_PHASE) &&
            (mPhaseTimer >= GetMDuration() + GetG1Duration() + GetSDuration() + GetG2Duration()) &&
            (mCurrentCellCyclePhase == M_PHASE)) 
        {
            mReadyToDivide = true;
        }
    }

    // If the cell is ready to divide, set the growth inhibited flag

    if (mpCell->GetCellData()->GetItem("FreeAreaFraction") < mpCell->GetCellData()->GetItem("p_beta") || mpCell->GetCellData()->GetItem("FreeSurfaceFraction") < mpCell->GetCellData()->GetItem("p_gamma"))
    {
        mpCell->GetCellData()->SetItem("growth inhibited", 1.0);
        mReadyToDivide = false;
    }
    else
    {
        mpCell->GetCellData()->SetItem("growth inhibited", 0.0);
    }

    if (mReadyToDivide && (mpCell->GetCellData()->GetItem("growth inhibited") == 0.0))
    {
        mpCell->GetCellData()->SetItem("Radius", (1/(2*std::sqrt(2))));
    }
    return mReadyToDivide;
}

void FixedDurationCellCycleModelWithContactInhibition::UpdateCellCyclePhase()
{
    double timeSinceBirth = GetAge();
    assert(timeSinceBirth >= 0);
    
    // Update cell growth phase timer
    // double change_in_cell_age = timeSinceBirth - mLastCellAge;
    // mPhaseTimer += change_in_cell_age;
    // mLastCellAge = timeSinceBirth;
    mPhaseTimer = mpCell->GetCellData()->GetItem("cell age");
    if (timeSinceBirth == 0.0)
    {
        mPhaseTimer = 0.0;
    }
    
    if (mPhaseTimer <= 0.0)
    {
        mPhaseTimer = timeSinceBirth;
    }
    else 
    {
        mPhaseTimer += (SimulationTime::Instance()->GetTimeStep())*(1.0 - mpCell->GetCellData()->GetItem("growth inhibited"));
        mpCell->SetBirthTime(SimulationTime::Instance()->GetTime() - mPhaseTimer);
    }
    // PRINT_3_VARIABLES(mPhaseTimer,mpCell->GetCellId(), mpCell->GetCellData()->GetItem("growth inhibited"));
    // PRINT_2_VARIABLES(mpCell->GetCellData()->GetItem("FreeAreaFraction") , mpCell->GetCellData()->GetItem("FreeSurfaceFraction"));

    SetPhaseTimer(mPhaseTimer);
    mpCell->GetCellData()->SetItem("cell age", mPhaseTimer);

    // double cells_radius = 0.3 + 0.2*(mPhaseTimer/(GetG1Duration() + GetSDuration() + GetG2Duration() + GetMDuration()));
    double cells_radius = (1/(2*std::sqrt(2)))*std::sqrt(1 + (mPhaseTimer/(GetG1Duration() + GetSDuration() + GetG2Duration() + GetMDuration())));
    mpCell->GetCellData()->SetItem("Radius", cells_radius);

    // Select the correct phase
    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
    else if (mPhaseTimer < GetG1Duration())
    {
        // if (mpCell->GetCellData()->GetItem("growth inhibited") != 0.0) {
        //     //std::cout << "Cell inhibited\n";
        //     mPhaseTimer -= change_in_cell_age;
        //     return;
        // }
        mCurrentCellCyclePhase = G_ONE_PHASE;
    }
    else if (mPhaseTimer <  GetG1Duration() + GetSDuration())
    {
        mCurrentCellCyclePhase = S_PHASE;
    }
    else if (mPhaseTimer < GetG1Duration() + GetSDuration() + GetG2Duration())
    {
        // if (mpCell->GetCellData()->GetItem("growth inhibited") != 0.0) {
        //     //std::cout << "Cell inhibited\n";
        //     mPhaseTimer -= change_in_cell_age;
        //     return;
        // }
        mCurrentCellCyclePhase = G_TWO_PHASE;
    }
    else if (mPhaseTimer < GetG1Duration() + GetSDuration() + GetG2Duration() + GetMDuration())
    {
        mCurrentCellCyclePhase = M_PHASE;
    }
    // PRINT_VARIABLE(mCurrentCellCyclePhase);

}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FixedDurationCellCycleModelWithContactInhibition)
