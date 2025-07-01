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

#ifndef TEST02MONOLAYERGROWTH_HPP_
#define TEST02MONOLAYERGROWTH_HPP_

#include <cxxtest/TestSuite.h> 

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp" 

#include "SmartPointers.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "MutableMesh.hpp"
#include "NodesOnlyMesh.hpp"
#include "HoneycombVertexMeshGenerator.hpp"

#include "CellsGenerator.hpp"

#include "BernoulliTrialWithContactInhibitionCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "CellProliferativeTypesCountWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellAncestorWriter.hpp"

#include "OffLatticeSimulation.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"

#include "VolumeTrackingModifier.hpp"
#include "SimpleTargetAreaModifier.hpp"
// #include "VertexBoundaryRefinementModifier.hpp"
// #include "SmoothVertexEdgesModifier.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"
#include "Debug.hpp"

#include "CicularityCalcModifier.hpp"
#include "BoundaryCellWriter.hpp"


/*
 *  This is where you can set parameters to be used in all the simulations.
 */

static const double M_END_TIME = 1; //100
static const double M_DT_TIME = 0.005;
static const double M_SAMPLE_TIME = 10;

static const double M_CONTACT_INHIBITION_LEVEL = 0.8;
static const double M_STEM_CELL_DIVISION_PROBABILITY = 0.1;
static const double M_STEM_CELL_MINIMUM_DIVISION_AGE = 1.0;

static const double M_INITIAL_WIDTH = 10;
static const double M_INITIAL_LENGTH = 10;

static const double M_DOMAIN_WIDTH = 5.0;
static const double M_DOMAIN_X_MIN = -5.0;
static const double M_DOMAIN_X_MAX = 5.0;
static const double M_DOMAIN_Y_MIN = -5.0;
static const double M_DOMAIN_Y_MAX = 5.0;

static const std::string M_HEAD_FOLDER = "OpenVT/Test02MonolayerGrowth";

class Test02MonolayerGrowth: public AbstractCellBasedWithTimingsTestSuite
{
private:

    /**
    * Helper method. Iterate over all cells and 'round out' the domain by
    * killing those cells whose centres are located outside a given region.
    * 
    * @param rCellPopulation a cell population
    */
    void RoundOutCellPopulation(AbstractCellPopulation<2>& rCellPopulation)
    {
        if (bool(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&rCellPopulation)))
        {
            std::set<unsigned> location_indices;
            std::set<unsigned> ghost_node_indices = dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&rCellPopulation)->GetGhostNodeIndices();

            for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
                cell_iter != rCellPopulation.rGetCells().end();)
            {
                // Get the coordinates of this cell centre
                c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                unsigned location_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
 
                // if ((fabs(y-x)>M_DOMAINWIDTH) && (x<M_DOMAIN_X_MIN) && (x>M_DOMAIN_X_MAX) && (y<M_DOMAIN_Y_MIN) && (y>M_DOMAIN_Y_MAX))
                if ( norm_2(centre_of_cell) > M_DOMAIN_WIDTH)
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

                // if ((fabs(y-x)>M_DOMAINWIDTH) && (x<M_DOMAIN_X_MIN) && (x>M_DOMAIN_X_MAX) && (y<M_DOMAIN_Y_MIN) && (y>M_DOMAIN_Y_MAX))
                if ( norm_2(centre_of_cell) > M_DOMAIN_WIDTH)
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
            else 
            {
              NEVER_REACHED;
            }
        }
    }

    /*
     * This is a helper method to generate cells and is used in all simulations.
     */ 

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells, double equilibriumVolume, double quiescentVolumeFraction, 
            double stemCellDivisionProbability, double stemCellMinimumDivisionAge)
    {

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_transit_cell_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_stem_cell_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        for (unsigned i=0; i<num_cells; i++)
        {
            BernoulliTrialWithContactInhibitionCellCycleModel* p_model = new BernoulliTrialWithContactInhibitionCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetEquilibriumVolume(equilibriumVolume);
            p_model->SetQuiescentVolumeFraction(quiescentVolumeFraction);
            p_model->SetStemCellDivisionProbability(stemCellDivisionProbability);
            p_model->SetStemCellMinimumDivisionAge(stemCellMinimumDivisionAge);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_cell_type);
            double ave_stem_cell_cycle_duration = p_model->GetAverageStemCellCycleTime();
            double birth_time = - RandomNumberGenerator::Instance()->ranf() * ave_stem_cell_cycle_duration;
            p_cell->SetBirthTime(birth_time);

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", sqrt(3.0)/2.0);

            rCells.push_back(p_cell);
        }
    }

public:

   /* 
     * == OS ==
     *
     * Simulate a growing monolayer using the
     * Overlapping Spheres model.
     */
    void TestNodeBasedGrowingMonolayer()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Node";

        /* 
         * == Default cut-off ==
         */
         // Create a simple mesh
        HoneycombMeshGenerator generator(2.0*M_INITIAL_WIDTH, 3.0*M_INITIAL_LENGTH,0);
        boost::shared_ptr<MutableMesh<2,2> > p_generating_mesh = generator.GetMesh();

        p_generating_mesh->Translate(-M_INITIAL_WIDTH,-0.75*sqrt(3.0)*M_INITIAL_LENGTH);

        double cut_off_length = 1.5; //this is the default

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(), cells, M_PI*0.25,M_CONTACT_INHIBITION_LEVEL,
                M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5
        /* If want to make default node the same size as mesh, need to set the mature volume to around
         * 0.92, r\approx 0.52.
         */

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();

        //cell_population.AddCellWriter<BoundaryCellWriter>();

        RoundOutCellPopulation(cell_population);

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_TIME); //50
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

        // TODO FIX THE CIRCULARITY CALCULATION
        // MAKE_PTR(CicularityCalcModifier<2>, circularity_modifier);
        // circularity_modifier->SetOutputDirectory(output_directory);
        // simulator.AddSimulationModifier(circularity_modifier);

        // Mark Ancestors
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation
        simulator.Solve();

        // Clear memory
        delete p_mesh;

    }

    /* 
     * == VT ==
     *
     * Simulate a growing monolayer using the
     * Bounded Voronoi Tesselation model.
     */
    void TestMeshBasedGrowingMonolayer()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Mesh";
         // Create mesh
        HoneycombMeshGenerator generator(2.0*M_INITIAL_WIDTH, 3.0*M_INITIAL_LENGTH);
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

        p_mesh->Translate(-M_INITIAL_WIDTH,-0.75*sqrt(3.0)*M_INITIAL_LENGTH);

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(), cells,sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: sqrt(3.0)/2.0
        // Note: Tissue shrinks without proliferation due to boundary effects

        // Create tissue
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        //cell_population.SetBoundVoronoiTessellation(true);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();

        //cell_population.AddCellWriter<BoundaryCellWriter>();

        RoundOutCellPopulation(cell_population);

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_TIME); //50
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        simulator.AddForce(p_linear_force);

        // TODO FIX THE CIRCULARITY CALCULATION
        // MAKE_PTR(CicularityCalcModifier<2>, circularity_modifier);
        // circularity_modifier->SetOutputDirectory(output_directory);
        // simulator.AddSimulationModifier(circularity_modifier);

        // Mark Ancestors
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation
        simulator.Solve();
    }

    /* 
     * == VM ==
     *
     * Simulate a growing monolayer using the
     * Cell Vertex model.
     */
    void TestVertexBasedGrowingMonolayer() 
    {
        std::string output_directory = M_HEAD_FOLDER + "/Vertex";

        // Create cells
        HoneycombVertexMeshGenerator mesh_generator(2.0*M_INITIAL_WIDTH, 3.0*M_INITIAL_LENGTH);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = mesh_generator.GetMesh();
        // p_mesh->SetCellRearrangementThreshold(0.1);

        p_mesh->Translate(-M_INITIAL_WIDTH,-0.75*sqrt(3.0)*M_INITIAL_LENGTH+0.25);

        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(), cells,sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); //mature_volume = 1.0

        // Create tissue
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();

        //cell_population.AddCellWriter<BoundaryCellWriter>();

        RoundOutCellPopulation(cell_population);

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_TIME); //50
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

        // Mark Ancestors
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // TODO FIX THE CIRCULARITY CALCULATION
        // MAKE_PTR(CicularityCalcModifier<2>, circularity_modifier);
        // circularity_modifier->SetOutputDirectory(output_directory);
        // simulator.AddSimulationModifier(circularity_modifier);

        // Add target area modifier
        // MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        // p_growth_modifier->SetGrowthDuration(0);
        // simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();
    }
};

#endif /* TEST02MONOLAYERGROWTH_HPP_ */
