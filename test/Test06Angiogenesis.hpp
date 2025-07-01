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
#include "DifferentialAdhesionPottsUpdateRule.hpp"
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
        
        double M_DOMAIN_WIDTH = 50.0;
        double M_DUDT_COEFICIENT = 1.0;
        double M_DIFFUSION_COEFICIENT = 1e2; //1e-3
        double M_CONSTANT_CELL_SECRETION = 1e-2; //1e-3
        double M_LINEAR_SECRETION = 1e-2; //1e-3
        double cell_area = 25.0;
        double box_h = 1.0;
        double dt = 0.01;
        double end_time = 1.0;
        double temperature = 1;


        EXIT_IF_PARALLEL;

        PottsMeshGenerator<2> generator((int) M_DOMAIN_WIDTH, 4, 4, (int) M_DOMAIN_WIDTH, 4, 4, 1, 1, 1, 4u);  // Parameters are: lattice sites across; num elements across; element width; lattice sites up; num elements up; element height; and element spacing.
        boost::shared_ptr<PottsMesh<2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
  
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetTemperature(temperature);
        cell_population.SetNumSweepsPerTimestep(1);
        
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("PottsBasedMonolayer");
        simulator.SetEndTime(end_time);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(10);

        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(25);
        
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        
        // MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_area_update_rule);
        // p_surface_area_update_rule->SetMatureCellTargetSurfaceArea(16);
        // p_surface_area_update_rule->SetDeformationEnergyParameter(0.5);
        // simulator.AddUpdateRule(p_surface_area_update_rule);

        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(40);
        p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(20);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        MAKE_PTR(VegfChemotaxisPottsUpdateRule<2>, p_chemotaxis_update_rule);
        //p_chemotaxis_update_rule->SetChemotaxisParameter(40);
        simulator.AddUpdateRule(p_chemotaxis_update_rule);

        // Set initial conditions
        cell_population.SetDataOnAllCells("VEGF", 0.0);

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-0.5,-0.5);
        ChastePoint<2> upper(M_DOMAIN_WIDTH-0.5, M_DOMAIN_WIDTH-0.5);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create PDE and BCS
        double dudt_coeficient = M_DUDT_COEFICIENT;
        double diffusion_coeficient = M_DIFFUSION_COEFICIENT;
        double constant_cell_uptake = M_CONSTANT_CELL_SECRETION*cell_area; 
        double linear_cell_uptake = M_LINEAR_SECRETION*cell_area;
        double constant_uptake = 0.0; 
        double linear_uptake = - M_LINEAR_SECRETION;
        // double constant_uptake = M_CONSTANT_SECRETION*cell_area; 
        // double linear_uptake = - M_LINEAR_CELL_SECRETION;
        bool is_volume_scaled = false;
        MAKE_PTR_ARGS(PottsParabolicPde<2>, p_pde, (cell_population, constant_cell_uptake, linear_cell_uptake, constant_uptake, linear_uptake, diffusion_coeficient, dudt_coeficient, is_volume_scaled));       
        //MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, constant_cell_uptake, linear_cell_uptake, constant_uptake, linear_uptake, diffusion_coeficient, dudt_coeficient, is_volume_scaled));       
        //MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, constant_cell_uptake, linear_cell_uptake, constant_uptake, linear_uptake, diffusion_coeficient, is_volume_scaled));       
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));


        // Create a PDE modifier and set the name of the dependent variable in the PDE
        bool is_neuman_bcs = false;
        
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
