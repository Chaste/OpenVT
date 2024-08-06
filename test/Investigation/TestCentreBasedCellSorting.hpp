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

#ifndef TESTCENTREBASEDCELLSORTING_HPP_
#define TESTCENTREBASEDCELLSORTING_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellLabel.hpp"
#include "MediumCellProperty.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "HeterotypicBoundaryLengthWriter.hpp"
#include "SortingCellTypesWriter.hpp"

#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "RandomMotionForce.hpp"
#include "SortingRandomMotionForce.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "ExtendedDifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "RandomMotionForce.hpp"

#include "OnLatticeSimulation.hpp"
#include "CellPopulationAdjacencyMatrixWriter.hpp"

#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "GhostNodeRemovalModifier.hpp"

#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

/*
 *  This is where you can set parameters to be used in all the simulations.
 *
 *  The first block (commented out) are the original parameter values.
 *  The second block are parameters for a much shorter simulation, and are used for continuous testing with Chaste.
 */

//static const double M_TIME_TO_STEADY_STATE = 10; //10
//static const double M_TIME_FOR_SIMULATION = 100; //100
//static const double M_NUM_CELLS_ACROSS = 20; //20 // this ^2 cells
//static const double M_CELL_FLUCTUATION = 1.0;

static const double M_TIME_TO_STEADY_STATE = 10.0; //10
static const double M_TIME_FOR_SIMULATION = 1010.0; //100
static const double M_NUM_CELLS_ACROSS = 10; //20 // this ^2 cells
static const double M_CELL_FLUCTUATION = 1.0;

static const std::string M_HEAD_FOLDER = "OpenVT/Prototyping";

class TestCentreBasedCellSorting : public AbstractCellBasedWithTimingsTestSuite
{
private:

    /*
     * This is a helper method to randomly label internal cells and make a boundary of 'GhostCells'.
     */ 

    void RandomlyLabelCells(AbstractCellPopulation<2>& rCellPopulation, boost::shared_ptr<AbstractCellProperty> pLabel, double labelledRatio, boost::shared_ptr<AbstractCellProperty> pMediumLabel)
    {
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
        {
            // Get the coordinates of this cell centre
            c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double x = centre_of_cell[0];
            double y = centre_of_cell[1];

            if ((x<0.5*M_NUM_CELLS_ACROSS) || (x>1.5*M_NUM_CELLS_ACROSS) || (y<0.866*0.5*M_NUM_CELLS_ACROSS) || (y>0.866*1.5*M_NUM_CELLS_ACROSS))
            {   
                (*cell_iter)->AddCellProperty(pMediumLabel);
            }
            else
            {
                if (RandomNumberGenerator::Instance()->ranf() < labelledRatio)
                {
                (*cell_iter)->AddCellProperty(pLabel);
                }
            }
        }
    }

public:

    /*
     * == OS ==
     *
     * Simulate a population of cells exhibiting cell sorting using the
     * Overlapping Sphere model.
     */
    void noTestNodeBasedMonolayerCellSorting()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Node";

        // Create a simple mesh
        HoneycombMeshGenerator generator(2*M_NUM_CELLS_ACROSS, 2*M_NUM_CELLS_ACROSS, 0);
        boost::shared_ptr<TetrahedralMesh<2,2> > p_generating_mesh = generator.GetMesh();

    		//Extended to allow sorting for longer distances
        double cut_off_length = 2.5; 

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, cut_off_length);

        // Set up cells, one for each Node
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<SortingCellTypesWriter>();
        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();
        cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);

        // Set time step and end time for simulation
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE);

        // Create a force law and pass it to the simulation
        MAKE_PTR(ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
        p_differential_adhesion_force->SetMeinekeSpringStiffness(50.0);
        p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1);
        p_differential_adhesion_force->SetLabelledMediumSpringConstantMultiplier(0.05);
        p_differential_adhesion_force->SetUnLabelledMediumSpringConstantMultiplier(0.1);
        p_differential_adhesion_force->SetMediumMediumSpringConstantMultiplier(0.1);

        p_differential_adhesion_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_differential_adhesion_force);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.05); //0.1 causes dissasociation, 0.001 is not enough
        simulator.AddForce(p_random_force);

        // Run simulation
        simulator.Solve();

        // Now label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        boost::shared_ptr<AbstractCellProperty> p_medium_state(CellPropertyRegistry::Instance()->Get<MediumCellProperty>());
        RandomlyLabelCells(simulator.rGetCellPopulation(), p_state, 0.5, p_medium_state);

        // Adjust parameters
        p_random_force->SetMovementParameter(0.05*M_CELL_FLUCTUATION); //0.1 causes dissasociation


        // Run simulation
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), M_NUM_CELLS_ACROSS*M_NUM_CELLS_ACROSS);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }
    /*
     * == VT ==
     *
     * Simulate a population of cells exhibiting cell sorting using the
     * Voronoi tesselation model.
     */
    void TestMeshBasedWithGhostsMonolayerCellSorting()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Mesh";

        // Create a simple mesh
        unsigned num_ghosts = 20;
        HoneycombMeshGenerator generator(2*M_NUM_CELLS_ACROSS, 2*M_NUM_CELLS_ACROSS, num_ghosts);
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

        // Set up cells, one for each non ghost Node
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_differentiated_type);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<SortingCellTypesWriter>();
        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);

        // Set time step and end time for simulation
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE);

        // Create a force law and pass it to the simulation
        MAKE_PTR(ExtendedDifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
        p_differential_adhesion_force->SetMeinekeSpringStiffness(50.0);
        p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1);
        p_differential_adhesion_force->SetLabelledMediumSpringConstantMultiplier(0.01);
        p_differential_adhesion_force->SetUnLabelledMediumSpringConstantMultiplier(0.1);
        p_differential_adhesion_force->SetMediumMediumSpringConstantMultiplier(0.1);
        simulator.AddForce(p_differential_adhesion_force);

        // Add some noise to avoid local minimum
        MAKE_PTR(SortingRandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(1.0); //0.1
        simulator.AddForce(p_random_force);

        // Add Modifier to remove internal ghost nodes
        MAKE_PTR(GhostNodeRemovalModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);   

        // Run simulation
        simulator.Solve();

        // Now label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        boost::shared_ptr<AbstractCellProperty> p_medium_state(CellPropertyRegistry::Instance()->Get<MediumCellProperty>());
        RandomlyLabelCells(simulator.rGetCellPopulation(), p_state, 0.5, p_medium_state);

        // Adjust parameters
        p_random_force->SetMovementParameter(0.1*M_CELL_FLUCTUATION); //0.1 causes dissasociation

        // Run simulation
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), M_NUM_CELLS_ACROSS*M_NUM_CELLS_ACROSS);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }
};

#endif /* TEST05CELLSORTING_HPP_ */
