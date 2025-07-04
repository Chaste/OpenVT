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

#ifndef TEST03WOUNDHEALING_HPP_
#define TEST03WOUNDHEALING_HPP_

#include <cxxtest/TestSuite.h> 

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp" 
#include "CheckpointArchiveTypes.hpp"

#include "SmartPointers.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "PeriodicNodesOnlyMesh.hpp"
#include "MutableMesh.hpp"
#include "Toroidal2dMesh.hpp"
#include "ToroidalHoneycombMeshGenerator.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"

#include "CellsGenerator.hpp"

#include "NoCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "CellVolumesWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "BoundaryNodeWriter.hpp"

#include "OffLatticeSimulation.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"

#include "SimpleTargetAreaModifier.hpp"
#include "VolumeTrackingModifier.hpp"

#include "PlaneBoundaryCondition.hpp"

//#include "GhostNodeRemovalModifier.hpp"
#include "VoidAreaModifier.hpp"
//#include "VertexBoundaryRefinementModifier.hpp"
//#include "CurvedVertexEdgesModifier.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp" 
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"
#include "Debug.hpp"

#include "BoundaryCellWriter.hpp"

// #include "CurvedVertexEdgesModifier.hpp"
// #include "SmoothVertexEdgesModifier.hpp"
// #include "JaggedVertexEdgesModifier.hpp"

/*
 *  This is where you can set parameters to be used in all the simulations.
 */

static const double M_END_STEADY_STATE = 0.01; //0.5;
static const double M_END_TIME = 0.02; //5;
static const double M_DT_TIME = 0.01;
static const double M_SAMPLE_TIME = 1; //10;

// Both Width and Length must be EVEN numbers here
static const double M_DOMAIN_WIDTH = 12;
static const double M_DOMAIN_LENGTH = 14;
static const double M_DOMAIN_SCALING = 0.8;
static const double M_PERIODIC_WIDTH = 10; //M_DOMAIN_WIDTH*M_DOMAIN_SCALING;
static const double M_PERIODIC_HEIGHT = 10; //M_DOMAIN_LENGTH*0.5*sqrt(3)*M_DOMAIN_SCALING;

static const double M_HOLEWIDTH = 2.0;
static const double M_HOLE_X_MIN = 2.0;
static const double M_HOLE_X_MAX = 8.0;
static const double M_HOLE_Y_MIN = 2.0;
static const double M_HOLE_Y_MAX = 8.0;

static const std::string M_HEAD_FOLDER = "OpenVT/Test03WoundHealing";



class Test03WoundHealing: public AbstractCellBasedWithTimingsTestSuite
{
private:
    /**
    * Helper method. Smooth out edges of a vertex mesh.
    * 
    * @param rCellPopulation a cell population
    */
    void SmoothVertexMeshEdges(AbstractCellPopulation<2>& rCellPopulation)
    {
        MutableVertexMesh<2, 2>& r_mesh = static_cast<VertexBasedCellPopulation<2>* >(&rCellPopulation)->rGetMesh();

        for (VertexMesh<2,2>::NodeIterator node_iter = r_mesh.GetNodeIteratorBegin();
            node_iter != r_mesh.GetNodeIteratorEnd();
            ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
            if (containing_element_indices.size() == 1)
            {
                // Get this element
                unsigned elem_index = (*containing_element_indices.begin());

                VertexElement<2,2>* p_element = r_mesh.GetElement(elem_index);

                // Remove node from this element and delete the node
                p_element->DeleteNode(p_element->GetNodeLocalIndex(node_index));
                r_mesh.DeleteNodePriorToReMesh(node_index);
            }
        }
        r_mesh.ReMesh();
    }

    /**
    * Helper method. Iterate over all cells and define the 'hole' by
    * killing those cells whose centres are located in a given region.
    * 
    * @param rCellPopulation a cell population
    */
    void CreateHoleInCellPopulation(AbstractCellPopulation<2>& rCellPopulation)
    {
        if (bool(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&rCellPopulation)))
        {
            std::set<unsigned> location_indices;
            std::set<unsigned> ghost_node_indices;

            for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
                cell_iter != rCellPopulation.rGetCells().end();)
            {
                // Get the coordinates of this cell centre
                c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                double x = centre_of_cell[0];
                double y = centre_of_cell[1];
                unsigned location_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
 
                if ((fabs(y-x)<M_HOLEWIDTH) && (x>M_HOLE_X_MIN) && (x<M_HOLE_X_MAX) && (y>M_HOLE_Y_MIN) && (y<M_HOLE_Y_MAX))
                {   
                    // Delete cell and store it as a ghost node
                    rCellPopulation.RemoveCellUsingLocationIndex(location_index, (*cell_iter));
                    // Update vector of cells
                    cell_iter = rCellPopulation.rGetCells().erase(cell_iter);
 
                    // Change to chost node            
                    ghost_node_indices.insert(location_index);
                }
                else
                {
                    ++cell_iter;
                }
            }
            dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&rCellPopulation)->SetGhostNodes(ghost_node_indices);
            //dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&rCellPopulation)->RemoveDeadCells();
            //dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&rCellPopulation)->Update();
        }
        else
        {
            for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                    cell_iter != rCellPopulation.End();
                    ++cell_iter)
            {
                // Get the coordinates of this cell centre
                c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                double x = centre_of_cell[0];
                double y = centre_of_cell[1];

                if ((fabs(y-x)<M_HOLEWIDTH) && (x>M_HOLE_X_MIN) && (x<M_HOLE_X_MAX) && (y>M_HOLE_Y_MIN) && (y<M_HOLE_Y_MAX))
                {   
                    cell_iter->Kill();
                }
            }
            
            /* Need to remove the deleted cells and call update note this is usually
            * performed in the Solve() method of the simulation class.
            */
            if (bool(dynamic_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation)))
            {
                dynamic_cast<NodeBasedCellPopulation<2>* >(&rCellPopulation)->RemoveDeadCells();
                dynamic_cast<NodeBasedCellPopulation<2>* >(&rCellPopulation)->Update();
            }
            else if (bool(dynamic_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation)))
            {
                dynamic_cast<MeshBasedCellPopulation<2>* >(&rCellPopulation)->RemoveDeadCells();
                dynamic_cast<MeshBasedCellPopulation<2>* >(&rCellPopulation)->Update();
            }
            else if (bool(dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation)))
            {
                dynamic_cast<VertexBasedCellPopulation<2>* >(&rCellPopulation)->RemoveDeadCells();
                dynamic_cast<VertexBasedCellPopulation<2>* >(&rCellPopulation)->Update();
            }
        }
    }

public:

    /* 
     * == OS ==
     *
     * Simulate an internal void using the
     * Overlapping Spheres model.
     *
     * Default Cut-off = 1.5
     */
    void TestNodeBasedInternalVoid()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Node/Pre-void";
        /* 
         * == Pre-void == 
         */
         // Create simple mesh
        HoneycombMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH, 0);
        boost::shared_ptr<MutableMesh<2, 2> > p_generating_mesh = generator.GetMesh();
        p_generating_mesh->Scale(M_DOMAIN_SCALING, M_DOMAIN_SCALING);

        double cut_off_length = 1.5; //this is the default

        // Convert this to a PeriodicNodesOnlyMesh
        c_vector<double,2> periodic_width = zero_vector<double>(2);
        periodic_width[0] = M_PERIODIC_WIDTH;
        periodic_width[1] = M_PERIODIC_HEIGHT;
        PeriodicNodesOnlyMesh<2>* p_mesh = new PeriodicNodesOnlyMesh<2>(periodic_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<BoundaryCellWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        p_linear_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_linear_force);

        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
        voidarea_modifier->SetOutputDirectory(output_directory);
        // voidarea_modifier->SetCutoff(cut_off_length);
        simulator.AddSimulationModifier(voidarea_modifier);

        // Run simulation
        simulator.Solve();

        // Save simulation in steady state
		    CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        // Clear memory
        delete p_mesh;

        // Load steady state
        OffLatticeSimulation<2>* p_simulator_1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
        NodeBasedCellPopulation<2>* p_cell_population_1 = static_cast<NodeBasedCellPopulation<2>*>(&(p_simulator_1->rGetCellPopulation()));

        std::string output_directory_1 =  M_HEAD_FOLDER + "/Node/Post-Void";

        // Now remove cells in a given region using a helper method
        CreateHoleInCellPopulation(*p_cell_population_1);
    
        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier1);
        voidarea_modifier1->SetOutputDirectory(output_directory_1);
        voidarea_modifier1->SetPixelSeparation(0.02);
        // voidarea_modifier1->SetCutoff(cut_off_length);
        p_simulator_1->AddSimulationModifier(voidarea_modifier1);
    
        // Reset timestep, sampling timestep and end time for simulation and run for a further duration
        p_simulator_1->SetDt(M_DT_TIME);
        p_simulator_1->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        p_simulator_1->SetEndTime(M_END_TIME);
        p_simulator_1->SetOutputDirectory(output_directory_1);
        p_simulator_1->Solve();

        // Tidy up
        delete p_simulator_1;

    }
   
    /* 
     * == VT ==
     * 
     * Simulate internal voide using the
     * Voronoi Tesselation model.
     */

    void TestMeshBasedInternalVoid()
    {
        std::string output_directory =  M_HEAD_FOLDER + "/Mesh/NoGhosts/Pre-Void";

        // Create mesh
        ToroidalHoneycombMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH, M_DOMAIN_SCALING, M_DOMAIN_SCALING);
        boost::shared_ptr<Toroidal2dMesh> p_mesh = generator.GetToroidalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create tissue
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<BoundaryCellWriter>();

        // Output Voroni for visualisation
        cell_population.AddPopulationWriter<VoronoiDataWriter>();
        cell_population.SetWriteVtkAsPoints(true);

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking Modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        p_linear_force->SetCutOffLength(1.0);
        simulator.AddForce(p_linear_force);

        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
        voidarea_modifier->SetOutputDirectory(output_directory);
        simulator.AddSimulationModifier(voidarea_modifier);

        // Run simulation
        simulator.Solve();

        // Save simulation in steady state
    		CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        // Load steady state
        OffLatticeSimulation<2>* p_simulator_2 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
        MeshBasedCellPopulation<2>* p_cell_population_2 = static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator_2->rGetCellPopulation()));

        std::string output_directory_2 =  M_HEAD_FOLDER + "/Mesh/Post-Void";

        // Now remove cells in a given region using a helper method
        CreateHoleInCellPopulation(*p_cell_population_2);
        
        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier_2);
        voidarea_modifier_2->SetOutputDirectory(output_directory_2);
        voidarea_modifier_2->SetPlotPixelContour(false);
        p_simulator_2->AddSimulationModifier(voidarea_modifier_2);
        
        // Bound the VT
        p_cell_population_2->SetBoundVoronoiTessellation(true);


        // Reset end time for simulation and run for a further duration
        p_simulator_2->SetOutputDirectory(output_directory_2);
        p_simulator_2->SetEndTime(M_END_TIME);
        p_simulator_2->Solve();

        // Tidy up
        delete p_simulator_2;
    }


    /* 
     * == VM ==
     * 
     * Simulation internal void using the
     * Cell Vertex model.
     */
    void TestVertexBasedInternalVoid()
    {
        std::string output_directory =  M_HEAD_FOLDER + "/Vertex/Pre-Void";

        /* 
         * == Pre-void == 
         */
         // Create mesh
        ToroidalHoneycombVertexMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH);
        boost::shared_ptr<Toroidal2dVertexMesh> p_mesh = generator.GetToroidalMesh();
        p_mesh->Scale(M_DOMAIN_SCALING, M_DOMAIN_SCALING);
        p_mesh->SetHeight(M_PERIODIC_HEIGHT);
        p_mesh->SetWidth(M_PERIODIC_WIDTH);
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Create tissue
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<BoundaryCellWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create Forces and pass to simulation NOTE : these are not the default ones and chosen to give a stable growing monolayer
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(2.0);
        simulator.AddForce(p_force);

        // Add target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetGrowthDuration(0);
        simulator.AddSimulationModifier(p_growth_modifier);
        
        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
        voidarea_modifier->SetOutputDirectory(output_directory);
        simulator.AddSimulationModifier(voidarea_modifier);
        
        // Smooth out edges to get nice box domain
        SmoothVertexMeshEdges(cell_population);
        
        // Run simulation
        simulator.Solve();

        // Save simulation in steady state
		    CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);
        
        // Load steady state
        OffLatticeSimulation<2>* p_simulator_2 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
        VertexBasedCellPopulation<2>* p_cell_population_2 = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator_2->rGetCellPopulation()));
        
        std::string output_directory_2 =  M_HEAD_FOLDER + "/Vertex/Post-Void";

        // Now remove cells in a given region using a helper method
        CreateHoleInCellPopulation(*p_cell_population_2);

        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier_2);
        voidarea_modifier_2->SetOutputDirectory(output_directory_2);
        p_simulator_2->AddSimulationModifier(voidarea_modifier_2);

        // MAKE_PTR(JaggedVertexEdgesModifier<2>, jagged_edge_modifier);
        // p_simulator_2->AddSimulationModifier(jagged_edge_modifier);

        // Reset timestep, sampling timestep and end time for simulation and run for a further duration
        p_simulator_2->SetDt(M_DT_TIME);
        p_simulator_2->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        p_simulator_2->SetEndTime(M_END_TIME);
        p_simulator_2->SetOutputDirectory(output_directory_2);
        p_simulator_2->Solve();

        // Tidy up
        delete p_simulator_2;
    
    }
};

#endif /* TEST03WOUNDHEALING_HPP_ */
