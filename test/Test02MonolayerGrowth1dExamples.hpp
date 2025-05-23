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

#ifndef TEST02MONOLAYERGROWTH1DEXAMPLES_HPP_
#define TEST02MONOLAYERGROWTH1DEXAMPLES_HPP_

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

#include "VolumeTrackingModifier.hpp"

#include "FixedDurationCellCycleModel.hpp"
#include "CellDataItemWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "TissueWidthWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "RandomCellKiller.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "RepulsionForce.hpp"

#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "GrowthInhibitionModifier.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "PottsMeshGeneratorExtended.hpp"
#include "OnLatticeSimulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "CellCentreLocationWriter.hpp"


#include "RandomNumberGenerator.hpp"

#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"


class Test02aMonlayerGrowth1dExamples : public AbstractCellBasedWithTimingsTestSuite
{
private:

    
public:

    void Test1dNodeChainCompression()
    {
        double end_time = 5.0;
        double dt = 0.01;
        unsigned output_timesteps = 1;
		double spring_stiffness = 50.0/3.0;
        //double compression = 0.5;
        //unsigned num_cells = 10;
        std::string base_type = "Test02MonlayerGrowthGrowth1d/Mesh";

        std::string tissue_types[2] = {"HomogeneousChain","HeterogeneousChain"};

        for (unsigned tissue_type_index = 0; tissue_type_index != 2; tissue_type_index++)
        {
            std::string tissue_type = tissue_types[tissue_type_index];

            std::string output_dir = base_type + "/" + tissue_type;

            std::string mesh_string;

            if(tissue_type.compare("HomogeneousChain")==0)
            {
                mesh_string = "projects/OpenVT/src/Test02MonolayerGrowth/1D_11_nodes";
            }
            else if(tissue_type.compare("HeterogeneousChain")==0)
            {
                mesh_string = "projects/OpenVT/src/Test02MonolayerGrowth/1D_11_plus_10_nodes";
            }
            else
            {
                NEVER_REACHED;
            }

            TrianglesMeshReader<1,1> mesh_reader_1d(mesh_string);
            MutableMesh<1,1> mesh;
            mesh.ConstructFromMeshReader(mesh_reader_1d);
            //mesh.Scale(compression);

            // Create cells
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            CellsGenerator<FixedDurationCellCycleModel, 1> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            MeshBasedCellPopulation<1> cell_population(mesh, cells);
            cell_population.AddCellWriter<CellCentreLocationWriter>();

            OffLatticeSimulation<1> simulator(cell_population);
            simulator.SetOutputDirectory(output_dir);
            simulator.SetDt(dt);
            simulator.SetSamplingTimestepMultiple(output_timesteps);
            simulator.SetEndTime(end_time);

            // Create a force law and pass it to the simulation
            MAKE_PTR(GeneralisedLinearSpringForce<1>, p_linear_force);
            p_linear_force->SetMeinekeSpringStiffness(spring_stiffness);
            p_linear_force->SetCutOffLength(1.5);
            simulator.AddForce(p_linear_force);

            // // Create a boundary condition for the left end 
            // c_vector<double,1> point = zero_vector<double>(1);
            // c_vector<double,1> normal = zero_vector<double>(1);
            // point(0) = 0.5;
            // normal(0) = -1.0;
            // MAKE_PTR_ARGS(PlaneBoundaryCondition<1>, p_bcs, (&cell_population, point, normal)); // y>0.5
            // simulator.AddCellPopulationBoundaryCondition(p_bcs);
        
            simulator.Solve();

            // Reset for next simulation
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }

    }




    void noTest1dPottsChainCompression()
    {
        unsigned start_index = 0;
        unsigned num_runs = 10;
        
        double end_time = 5.0;//1.0;
		double dt = 0.025; // EDITED TO MAKE HIT 9CD AT 1HR
        unsigned output_timesteps = 4;
        
        double target_area = 50.0;
        double target_area_parameter = 20.0;
        
        double temperature = 20;

        unsigned domain_length = 250;
        unsigned domain_width = 5;
        unsigned initial_cell_length = 5;

        //double compression = 0.5;
        unsigned num_cells = 11;
        unsigned num_edge_cells = 5;
       
        std::string base_type = "Test02MonlayerGrowthGrowth1d/Potts";

        std::string tissue_types[2] = {"HomogeneousChain","HeterogeneousChain"};

         // Loop over the random seed.
		for(unsigned sim_index=start_index; sim_index < start_index + num_runs; sim_index++)
		{
			std::cout << " Run number " << sim_index << "... \n" << std::flush;

			// Reseed the random number generator
			RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
			p_gen->Reseed(sim_index);

			std::stringstream sim_index_ss;
			sim_index_ss << "Run" << sim_index;
            std::string sim_index_string = sim_index_ss.str();
        
            for (unsigned tissue_type_index = 0; tissue_type_index != 2; tissue_type_index++)
            {
                std::string tissue_type = tissue_types[tissue_type_index];

                std::string output_dir = base_type + "/" + sim_index_string + "/" + tissue_type;

                std::string mesh_string;

                boost::shared_ptr<PottsMesh<2> > p_mesh;

                if(tissue_type.compare("HomogeneousChain")==0)
                {
                    std::vector<unsigned> cell_lengths(num_cells, initial_cell_length);

                    PottsMeshGeneratorExtended<2> generator(domain_length, num_cells, cell_lengths, domain_width, 1, domain_width,  1, 1, 1, 0, false, false, true, false);  
                    p_mesh = generator.GetMesh();
                }
                else if(tissue_type.compare("HeterogeneousChain")==0)
                {
                    std::vector<unsigned> cell_lengths;
                    for (unsigned i=0; i<num_edge_cells; i++)
                    {
                        cell_lengths.push_back(2.0*initial_cell_length);
                    }
                    for (unsigned i=0; i<num_cells; i++)
                    {
                        cell_lengths.push_back(initial_cell_length);
                    }
                    for (unsigned i=0; i<num_edge_cells; i++)
                    {
                        cell_lengths.push_back(2.0*initial_cell_length);
                    }

                    PottsMeshGeneratorExtended<2> generator(domain_length, num_cells + 2*num_edge_cells, cell_lengths, domain_width, 1, domain_width,  1, 1, 1, 0, false, false, true, false);
                    p_mesh = generator.GetMesh();
                }
                else
                {
                    NEVER_REACHED;
                }

                // Create cells
                std::vector<CellPtr> cells;
                MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
                CellsGenerator<FixedDurationCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

                PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.SetTemperature(temperature);
                cell_population.SetNumSweepsPerTimestep(1);
                cell_population.AddCellWriter<CellCentreLocationWriter>();

                OnLatticeSimulation<2> simulator(cell_population);
                simulator.SetOutputDirectory(output_dir);
                simulator.SetEndTime(end_time);
                simulator.SetDt(dt);
                simulator.SetSamplingTimestepMultiple(output_timesteps);

                MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
                p_volume_constraint_update_rule->SetMatureCellTargetVolume(target_area);
                p_volume_constraint_update_rule->SetDeformationEnergyParameter(target_area_parameter);
                simulator.AddUpdateRule(p_volume_constraint_update_rule);
                
                // MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_area_update_rule);
                // p_surface_area_update_rule->SetMatureCellTargetSurfaceArea(M_TARGET_CELL_SURFACE_AREA);
                // p_surface_area_update_rule->SetDeformationEnergyParameter(0.5);
                // simulator.AddUpdateRule(p_surface_area_update_rule);

                MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
                p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(20);
                p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(20);
                simulator.AddUpdateRule(p_adhesion_update_rule);
            
                simulator.Solve();

                // Reset for next simulation
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }
        }

    }


};

#endif /* TEST02MONOLAYERGROWTH1DEXAMPLES_HPP_ */
