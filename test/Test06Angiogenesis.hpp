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

#ifndef TEST06ANGIOGENESIS_HPP_
#define TEST06ANGIOGENESIS_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"
#include "PottsMeshGenerator.hpp"
#include "OnLatticeSimulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "VegfChemotaxisPottsUpdateRule.hpp"
#include "ExtendedAdhesionPottsUpdateRule.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "CellLabelWriter.hpp"

// #include "AveragedSourceEllipticPde.hpp"
// #include "EllipticPottsPdeModifier.hpp"
#include "PottsParabolicPde.hpp"
#include "ParabolicPottsPdeModifier.hpp"

class Test06Angiogenesis: public AbstractCellBasedWithTimingsTestSuite
{
private:

public:

    /*
        Simple Potts Monolayer
     */
    void TestMonolayer()
    {
        
        double M_DOMAIN_WIDTH = 100.0; //200
        double M_DUDT_COEFICIENT = 1.0;
        double M_DIFFUSION_COEFICIENT = 1.0; //1e-4; //1e-3
        double M_CONSTANT_CELL_SECRETION = 0.3; //1e-2; //1e-3
        double M_LINEAR_SECRETION = 0.3; //1e0; //1e-3
        double M_TARGET_CELL_VOLUME = 50;
        double M_CELL_VOLUME_PARAMETER = 5;
        //double M_CHEMOTAXIS_PARAMETER = 500; //50000
        double M_CELLCELL_ADHESION = 2;
        double M_CELLMEMBRANE_ADHESION = 8;
        
        // double M_TARGET_CELL_SURFACE_AREA = 16;
        double dt = 1.0/10.0;
        double end_time = 10.0;
        unsigned output_timesteps = 10;
        double temperature = 5; //25
        unsigned num_sweeps = 1;

        EXIT_IF_PARALLEL;

        unsigned initial_cell_width = 5u;
        unsigned initial_separation = 0u;//5u; //2u
        unsigned initial_num_cells = 10u; //20 //16u; // 32u

        PottsMeshGenerator<2> generator((int) M_DOMAIN_WIDTH, initial_num_cells, initial_cell_width, (int) M_DOMAIN_WIDTH, initial_num_cells, initial_cell_width, 1, 1, 1, initial_separation, false, true, true, false);  
        boost::shared_ptr<PottsMesh<2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
  
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetTemperature(temperature);
        cell_population.SetNumSweepsPerTimestep(num_sweeps);
        
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("PottsBasedMonolayer");
        simulator.SetEndTime(end_time);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(output_timesteps);

        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(M_TARGET_CELL_VOLUME);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(M_CELL_VOLUME_PARAMETER);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        
        // MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_area_update_rule);
        // p_surface_area_update_rule->SetMatureCellTargetSurfaceArea(M_TARGET_CELL_SURFACE_AREA);
        // p_surface_area_update_rule->SetDeformationEnergyParameter(0.5);
        // simulator.AddUpdateRule(p_surface_area_update_rule);

        // MAKE_PTR(ExtendedAdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(M_CELLCELL_ADHESION);
        p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(M_CELLMEMBRANE_ADHESION);
        // p_adhesion_update_rule->SetCellBoundaryNodeAdhesionEnergyParameter(100);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // MAKE_PTR(VegfChemotaxisPottsUpdateRule<2>, p_chemotaxis_update_rule);
        // p_chemotaxis_update_rule->SetChemotaxisParameter(M_CHEMOTAXIS_PARAMETER);
        // simulator.AddUpdateRule(p_chemotaxis_update_rule);

        // Set initial conditions
        cell_population.SetDataOnAllCells("VEGF", 0.0);

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(0.0,0.0);
        ChastePoint<2> upper(M_DOMAIN_WIDTH-1, M_DOMAIN_WIDTH-1);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create PDE and BCS
        double dudt_coeficient = M_DUDT_COEFICIENT;
        double diffusion_coeficient = M_DIFFUSION_COEFICIENT;
        double constant_cell_uptake = M_CONSTANT_CELL_SECRETION; 
        double linear_cell_uptake = M_LINEAR_SECRETION;
        double constant_uptake = 0.0; 
        double linear_uptake = - M_LINEAR_SECRETION;
        // double constant_uptake = M_CONSTANT_SECRETION*cell_area; 
        // double linear_uptake = - M_LINEAR_CELL_SECRETION;
        // bool is_volume_scaled = false;
        MAKE_PTR_ARGS(PottsParabolicPde<2>, p_pde, (cell_population, constant_cell_uptake, linear_cell_uptake, constant_uptake, linear_uptake, diffusion_coeficient, dudt_coeficient));       
        //MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, constant_cell_uptake, linear_cell_uptake, constant_uptake, linear_uptake, diffusion_coeficient, dudt_coeficient, is_volume_scaled));       
        //MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, constant_cell_uptake, linear_cell_uptake, constant_uptake, linear_uptake, diffusion_coeficient, is_volume_scaled));       
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        bool is_neuman_bcs = true;
        double box_h = 1.0; // to match with PottsMesh
        MAKE_PTR_ARGS(ParabolicPottsPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, is_neuman_bcs, p_cuboid, box_h));
        //MAKE_PTR_ARGS(EllipticPottsPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, is_neuman_bcs, p_cuboid, box_h));
        p_pde_modifier->SetDependentVariableName("VEGF");
        p_pde_modifier->SetBcsOnBoxBoundary(true);
        p_pde_modifier->SetOutputSolutionAtPdeNodes(true); 
        //p_pde_modifier->SetSolutionInterval(10); // Speed simulations up a bit
        //p_pde_modifier->SetOutputGradient(true);

        simulator.AddSimulationModifier(p_pde_modifier); 


        simulator.Solve();
    }
};

#endif /* TEST06ANGIOGENESIS_HPP_ */
