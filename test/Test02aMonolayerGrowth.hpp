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

#ifndef TEST02AMONOLAYERGROWTH_HPP_
#define TEST02AMONOLAYERGROWTH_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "DefaultCellProliferativeType.hpp"

#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "VolumeTrackingModifier.hpp"

#include "FixedDurationCellCycleModelWithGrowthInhibition.hpp"
#include "CellDataItemWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "TissueWidthWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "RandomCellKiller.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForceWithMinDistanceItem.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "RepulsionForce.hpp"

#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "GrowthInhibitionModifier.hpp"

#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"

#include "CaBasedCellPopulation.hpp"
#include "DiffusionCaUpdateRule.hpp"

#include "RandomNumberGenerator.hpp"

#include "PetscSetupAndFinalize.hpp"


class Test02aMonlayerGrowth : public AbstractCellBasedWithTimingsTestSuite
{
private:

    /*
     * This is a helper method to generate cells and is used in all simulations.
     */ 
    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells, bool randomiseBirthTime)
    {
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        for (unsigned i=0; i<num_cells; i++)
        {
            //UniformlyDistributedCellCycleModel* p_cycle_model = new UniformlyDistributedCellCycleModel();
            FixedDurationCellCycleModelWithGrowthInhibition* p_cycle_model = new FixedDurationCellCycleModelWithGrowthInhibition();
            p_cycle_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);

            double birth_time = 0.0;
            if (randomiseBirthTime) {
                birth_time = -RandomNumberGenerator::Instance()->ranf() * 18.0;
            }
            p_cell->SetBirthTime(birth_time);
            p_cycle_model->SetPhaseTimer(birth_time);


            p_cell->InitialiseCellCycleModel();

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", 1.0);
            p_cell->GetCellData()->SetItem("growth inhibited", 0.0);
            p_cell->GetCellData()->SetItem("Radius", 0.1);
            p_cell->GetCellData()->SetItem("cell age", birth_time);
            rCells.push_back(p_cell);
        }
     }

public:

    /*
     * Simulate growth of a tissue monolayer without diffusion. Starts with a single cell
     */
    void Test2DMonolayerWithoutDiffusionSingleCell()
    {
        static const double end_time = 28*24; // 28 days first 14 days and second 14 days can be separated 


        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        std::vector<double> center = {0.0, 0.0};
        double cut_off_length = 1.5; //this is the default
        Node<2> node(0, center.data(), false);
        p_mesh->AddNode(&node);
        p_mesh->SetMaximumInteractionDistance(cut_off_length);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells,false);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddPopulationWriter<TissueWidthWriter>();
        cell_population.SetUseVariableRadii(true);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Test02aMonlayerGrowth");
        simulator.SetDt(0.05);
        simulator.SetSamplingTimestepMultiple(20); // Every hour
        simulator.SetEndTime(end_time);

        simulator.SetOutputDivisionLocations(true);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForceWithMinDistanceItem<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(2.70);
        p_linear_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_linear_force);
        
        MAKE_PTR(GrowthInhibitionModifier<2>, p_growth_inhibition_modifier);
        simulator.AddSimulationModifier(p_growth_inhibition_modifier);


        simulator.Solve();

    }

};

#endif /* TEST02MONOLAYERGROWTH_HPP_ */
