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

#include "FixedGrowthModelWithContactInhibition_mod.hpp"
#include "Debug.hpp"

#include <cmath>

#include <iostream>

#include "RandomNumberGenerator.hpp"

template<class Archive>
void FixedGrowthModelWithContactInhibition_mod::serialize(Archive & archive, const unsigned int version)
{
    // Archive cell-cycle model using serialization code from AbstractSimplePhaseBasedCellCycleModel
    archive & boost::serialization::base_object<AbstractSimplePhaseBasedCellCycleModel>(*this);
}

FixedGrowthModelWithContactInhibition_mod::FixedGrowthModelWithContactInhibition_mod()
    : AbstractSimplePhaseBasedCellCycleModel()
{
    // SetPhaseDurations();
    mPhaseTimer = 0.0;
    mLastCellAge = 0.0;
    // mBeta = 1.0;  // Default threshold for area fraction below which growth is inhibited
    // mGamma = 1.0; // Default threshold for surface area fraction below which growth is inhibited
}

void FixedGrowthModelWithContactInhibition_mod::ResetForDivision()
{
    AbstractCellCycleModel::ResetForDivision();
    mCurrentCellCyclePhase = G_ONE_PHASE;
    mPhaseTimer = 0.0;
    mLastCellAge = 0.0;
}

// void FixedGrowthModelWithContactInhibition_mod::SetPhaseDurations()
// {
//     SetStemCellG1Duration(2.0);  
//     SetTransitCellG1Duration(2.0);
//     SetSDuration(3.0);
//     SetG2Duration(3.0);
//     SetMDuration(2.0);
//     SetMinimumGapDuration(3.0);
// }

void FixedGrowthModelWithContactInhibition_mod::SetG1Duration()
{
    assert(mpCell != NULL);  // Make sure cell exists

    mG1Duration = 7.0;
}

AbstractCellCycleModel* FixedGrowthModelWithContactInhibition_mod::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    FixedGrowthModelWithContactInhibition_mod* pCellCycleModel = new FixedGrowthModelWithContactInhibition_mod();
    return pCellCycleModel;
}

void FixedGrowthModelWithContactInhibition_mod::SetPhaseTimer(const double value) {
    mPhaseTimer = value;
}

double FixedGrowthModelWithContactInhibition_mod::GetPhaseTimer()
{
    return mPhaseTimer;
}

bool FixedGrowthModelWithContactInhibition_mod::ReadyToDivide() 
{
    // TRACE("UPDATING CELL CYCLE PHASE");
    // PRINT_VARIABLE((SimulationTime::Instance()->GetTime()) );

    assert(mpCell != nullptr);
    if (!mReadyToDivide)
    {   
        UpdateCellCyclePhase();
        if ((mCurrentCellCyclePhase == M_PHASE)) 
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
        mpCell->GetCellData()->SetItem("Radius", (1.0/(std::sqrt(M_PI))));
        mpCell->GetCellData()->SetItem("Deformable Radius", (1.0/(std::sqrt(M_PI))));
        mpCell->GetCellData()->SetItem("cell age", 0.0);

    }
    // TRACE("FINISHED UPDATE CELL CYCLE PHASE");

    return mReadyToDivide;

}

void FixedGrowthModelWithContactInhibition_mod::UpdateCellCyclePhase()
{
    
    double timeSinceBirth = GetAge();
    assert(timeSinceBirth >= 0);
    
    mPhaseTimer = mpCell->GetCellData()->GetItem("cell age");

    // At birth, set target area and birth age based on normal random variable
    if (timeSinceBirth <= (SimulationTime::Instance()->GetTimeStep()))
    {
        double initial_radius = mpCell->GetCellData()->GetItem("Radius");
        // double deformable_radius = mpCell->GetCellData()->GetItem("Deformable Radius");
        // mpCell->GetCellData()->SetItem("Initial_Radius", initial_radius);
        double parent_cell_target_area = mpCell->GetCellData()->GetItem("target area");
        double growth_rate = mpCell->GetCellData()->GetItem("growth rate");
        mPhaseTimer = 0.0;

        // generate normal random variable with mean 2 and std dev 0.4^2
        double final_area_i = -1.0;
        while (final_area_i <= 0.0)
        {
            final_area_i = RandomNumberGenerator::Instance()->NormalRandomDeviate(2.0, 0.4*0.4);
        }
        double initial_area = M_PI*initial_radius*initial_radius;
        double birth_time = (final_area_i - initial_area)/growth_rate;

        // double birth_time = (final_area_i - 0.5*parent_cell_target_area)/growth_rate;
        mpCell->GetCellData()->SetItem("birth age", birth_time);
        mpCell->GetCellData()->SetItem("target area", final_area_i);
        mpCell->GetCellData()->SetItem("Deformable Radius", (1.0/(std::sqrt(M_PI))));
        
    }
    
    if (mPhaseTimer <= 0.0)
    {
        mPhaseTimer = timeSinceBirth;
    }
    else 
    {
        int is_growth_inhibited = int (mpCell->GetCellData()->GetItem("growth inhibited") == 1.0);
        int is_big_enough = int (mpCell->GetCellData()->GetItem("Deformable Radius")  >= 0.999*mpCell->GetCellData()->GetItem("Radius") );
        is_big_enough = 1.0;
        
        mPhaseTimer += (SimulationTime::Instance()->GetTimeStep())*(1.0 - is_growth_inhibited)*(is_big_enough);
        mpCell->SetBirthTime(SimulationTime::Instance()->GetTime() - mPhaseTimer);
    }

    double cell_birth_age = mpCell->GetCellData()->GetItem("birth age");
    double target_area = mpCell->GetCellData()->GetItem("target area");
    
    mPhaseTimer = std::min(mPhaseTimer, cell_birth_age);
    SetPhaseTimer(mPhaseTimer);
    mpCell->GetCellData()->SetItem("cell age", mPhaseTimer);
    // if (cell_birth_age >= )

    // double cells_radius = (1/(2*std::sqrt(2)))*std::sqrt(1 + (mPhaseTimer/(cell_birth_age)));
    // double cells_radius = (0.5*std::sqrt(target_area/M_PI))*std::sqrt(1 + (mPhaseTimer/(cell_birth_age)));
    double cells_radius = (std::sqrt(0.5*target_area/M_PI))*std::sqrt(1 + (mPhaseTimer/(cell_birth_age)));

    mpCell->GetCellData()->SetItem("Radius", cells_radius);
    

    // Select the correct phase
    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
    else if (mPhaseTimer >= cell_birth_age)
    {
        double cells_radius = (std::sqrt(0.5*target_area/M_PI))*std::sqrt(2);
        mpCell->GetCellData()->SetItem("Radius", cells_radius);

        if (mpCell->GetCellData()->GetItem("Deformable Radius")  >= 0.999*mpCell->GetCellData()->GetItem("Radius") )
        {
            mCurrentCellCyclePhase = M_PHASE;
        }
        else
        {
            mCurrentCellCyclePhase = G_ZERO_PHASE;
        }
        // mCurrentCellCyclePhase = M_PHASE;

    }
    else
    {
        mCurrentCellCyclePhase = G_ONE_PHASE;
    }

}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FixedGrowthModelWithContactInhibition_mod)
