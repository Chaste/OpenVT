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

#ifndef MORPHOGENDEPENDENTCELLCYCLEMODEL_HPP_
#define MORPHOGENDEPENDENTCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"

/**
 * Cell cycle model which implements the proliferation model from:
 *
 * Smith et al Incorporating chemical signalling factors into cell-based models of growing epithelial tissues.
 * doi: 10.1007/s00285-011-0464-y
 */

class MorphogenDependentCellCycleModel : public AbstractCellCycleModel
{
private:

    friend class boost::serialization::access;

    /**
     * Boost Serialization method for archiving/checkpointing
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        archive & mGrowthRate;
        archive & mCurrentMass;
        archive & mTargetMass;
        archive & mMorphogenInfluence;
    }

    /**
     * Growth rate of the cell, randomly assigned on birth
     */
    double mGrowthRate;

    /**
     * Stores the current mass of the cell, (Referred to as volume in the paper)
     */
    double mCurrentMass;

    // Parameters

    /**
     * Target mass of the cell (Referred to as volume in the paper)
     */
    double mTargetMass;

    /**
     * the level of influence the morphogen has on growth
     */
    double mMorphogenInfluence;

protected:

    /**
     * Protected copy-constructor for use by CreateCellCycleModel.
     * The only way for external code to create a copy of a cell cycle model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent cell cycle model will have had ResetForDivision() called just before CreateCellCycleModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    MorphogenDependentCellCycleModel(const MorphogenDependentCellCycleModel& rModel);

public:

    /**
     * Constructor.
     */
    MorphogenDependentCellCycleModel();

    /**
     * Overridden ReadyToDivideMethod
     *
     * @return whether the cell is ready to divide (enter M phase).
     *
     * Here we divide randomly basd on the size of the cell compared to the target size.
     *
     */
    virtual bool ReadyToDivide();

    /**
     * Overriden ResetForDivision method to halve the current cell size.
     *
     * Should only be called by the Cell Divide() method.
     */
    virtual void ResetForDivision();

    /**
     * Overridden builder method to create new instances of
     * the cell-cycle model.
     *
     * @return new cell-cycle model
     *
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * @param thresholdMorphogen
     */
    void SetCurrentMass(double currentMass);

    /**
     * @return mCurrentMass
     */
    double GetCurrentMass();

    /**
     * @param targetMass
     */
    void SetTargetMass(double targetMass);

    /**
     * @return mTargetMass
     */
    double GetTargetMass();

    /**
     * @param morphogenInfluence
     */
    void SetMorphogenInfluence(double morphogenInfluence);

    /**
     * @return mMorphogenInfluence
     */
    double GetMorphogenInfluence();

    /**
     * Helper method to generate the mGrowthRate from a truncated normal distribution.
     */
    void GenerateGrowthRate();

    /**
     * Overridden GetAverageTransitCellCycleTime() method.
     * @return time
     */
    double GetAverageTransitCellCycleTime();

    /**
     * Overridden GetAverageStemCellCycleTime() method.
     * @return time
     */
    double GetAverageStemCellCycleTime();

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MorphogenDependentCellCycleModel)

#endif // MORPHOGENDEPENDENTCELLCYCLEMODEL_HPP_
