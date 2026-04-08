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
#include "NodeLocationWriter.hpp"

#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "VolumeTrackingModifier.hpp"

#include "CellDataItemWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "TissueWidthWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "RandomCellKiller.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "RepulsionForce.hpp"

#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "FixedDurationCellCycleModelWithGrowthInhibition.hpp"
#include "FixedDurationCellCycleModelWithContactInhibition.hpp"
#include "GrowthInhibitionModifier.hpp"

#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"

#include "CaBasedCellPopulation.hpp"
#include "DiffusionCaUpdateRule.hpp"

#include "RandomNumberGenerator.hpp"
#include "OffLatticeSimulationWithPopulationBasedStoppingEvent.hpp"
#include "ForwardEulerNumericalMethod.hpp"

#include "PetscSetupAndFinalize.hpp"
// #include "CellContactInhibitionWriter.hpp"
#include "Debug.hpp"

#include "RandomMotionForce.hpp"


#include "FixedGrowthModelWithContactInhibition.hpp"
#include "ContactInhibitionModifier.hpp"
#include "GeneralisedLinearSpringForceWithMinDistanceItem.hpp"

// #include "FixedGrowthModelWithContactInhibition_mod.hpp"
// #include "ContactInhibitionModifier_mod.hpp"
// #include "GeneralisedLinearSpringForceWithMinDistanceItem_mod.hpp"


class Test02aMonlayerGrowth : public AbstractCellBasedWithTimingsTestSuite
{
private:

    /*
     * This is a helper method to generate cells and is used in all simulations.
     */ 
    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells, bool randomiseBirthTime, double p_beta, double p_gamma)
    {
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        for (unsigned i=0; i<num_cells; i++)
        {
            //UniformlyDistributedCellCycleModel* p_cycle_model = new UniformlyDistributedCellCycleModel();
            // FixedDurationCellCycleModelWithGrowthInhibition* p_cycle_model = new FixedDurationCellCycleModelWithGrowthInhibition();
            
            // FixedGrowthModelWithContactInhibition_mod* p_cycle_model = new FixedGrowthModelWithContactInhibition_mod();
            FixedGrowthModelWithContactInhibition* p_cycle_model = new FixedGrowthModelWithContactInhibition();
            p_cycle_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_stem_type);
      

            // generate normal random variable with mean 2 and std dev 0.4^2
            double growth_rate = 1.0/5.0;
            double final_area_i = -1.0;
            while (final_area_i <= 0.0)
            {
                final_area_i = RandomNumberGenerator::Instance()->NormalRandomDeviate(2.0, 0.4*0.4);
                // PRINT_VARIABLE(final_area_i);
            }
            double initial_area = 0.5*final_area_i;
            double birth_time = (final_area_i - initial_area)/growth_rate;


            p_cell->SetBirthTime(birth_time);
            // p_cycle_model->SetPhaseTimer(birth_time);
            // p_cycle_model->SetFreeAreaFraction(p_beta);
            // p_cycle_model->SetFreeSurfaceFraction(p_gamma);

            p_cell->InitialiseCellCycleModel();

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", final_area_i);
            // p_cell->GetCellData()->SetItem("growth rate", 0.2); //0.1
            p_cell->GetCellData()->SetItem("growth rate", growth_rate); //0.9116

            p_cell->GetCellData()->SetItem("birth age", 10.0);

            p_cell->GetCellData()->SetItem("growth inhibited", 0.0);
            p_cell->GetCellData()->SetItem("Radius", initial_area);
            p_cell->GetCellData()->SetItem("Deformable Radius", initial_area);
            // p_cell->GetCellData()->SetItem("Initial_Radius", 0.35);
            p_cell->GetCellData()->SetItem("cell age", birth_time);
            p_cell->GetCellData()->SetItem("FreeSurfaceFraction", 1.0);
            p_cell->GetCellData()->SetItem("FreeAreaFraction", 1.0);
            p_cell->GetCellData()->SetItem("p_beta", p_beta);
            p_cell->GetCellData()->SetItem("p_gamma", p_gamma);
            rCells.push_back(p_cell);
        }
     }

public:

    /*
     * Simulate growth of a tissue monolayer without diffusion. Starts with a single cell
     */
    void Test2DMonolayerWithoutDiffusionSingleCell()
    {

        double spring_stiffness = 0.0;
        // create a string identifier to take the string values of either linear, quadratic or log force law
        // std::string force_law = "quadratic"; // "log"; // "quadratic"; // "linear";
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-p_force_law"));
        std::string force_law = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-p_force_law");
        if (force_law == "linear")
        {
            spring_stiffness = 18.2816648;
        }
        else if (force_law == "quadratic")
        {
            // spring_stiffness = 88.54;
            spring_stiffness = 155.0;
        }
        else if (force_law == "log")
        {
            spring_stiffness = 15.6556;
        }
        PRINT_VARIABLE(force_law);
        PRINT_VARIABLE(spring_stiffness);

        // reset the random number generator
        double random_seed;
        try
        {
          TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-random_seed"));
          random_seed = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-random_seed");
        }
        catch(const std::exception& e)
        {
          random_seed = 0.0;
        }
        RandomNumberGenerator::Instance()->Reseed(random_seed);

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-p_beta"));
        double p_beta = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-p_beta");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-p_gamma"));
        double p_gamma = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-p_gamma");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-end_time"));
        double end_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-end_time");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-sample_rate"));
        double sample_rate = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-sample_rate");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-output_name"));
        std::string output_name = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-output_name");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-cut_off_length"));
        double cut_off_length = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-cut_off_length");


        // static const double end_time = 200; //28*24; // 28 days first 14 days and second 14 days can be separated 


        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        
        std::vector<double> cell_1 = {0.0, 0.0};
        Node<2> node1(0, cell_1.data(), false);
        p_mesh->AddNode(&node1);

        // std::vector<double> cell_2 = {0.5, 0.0};
        // Node<2> node2(1, cell_2.data(), false);
        // p_mesh->AddNode(&node2);

        // std::vector<double> cell_3 = {0.0, 0.5};
        // Node<2> node3(2, cell_3.data(), false);
        // p_mesh->AddNode(&node3);
      
        // std::vector<double> cell_4 = {-0.5, 0.0};
        // Node<2> node4(3, cell_4.data(), false);
        // p_mesh->AddNode(&node4);

        // std::vector<double> cell_5 = {0.0, -0.5};
        // Node<2> node5(4, cell_5.data(), false);
        // p_mesh->AddNode(&node5);

        p_mesh->SetMaximumInteractionDistance(1.1*cut_off_length);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells, false, p_beta, p_gamma);

        double division_separation = 2.0*(std::sqrt(1.0/M_PI));

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        // cell_population.AddCellWriter<CellContactInhibitionWriter>();
        cell_population.AddPopulationWriter<TissueWidthWriter>();
        cell_population.SetMeinekeDivisionSeparation(division_separation);
        cell_population.SetUseVariableRadii(true);

        OffLatticeSimulationWithPopulationBasedStoppingEvent simulator(cell_population);
        simulator.SetOutputDirectory(output_name);
        simulator.SetDt(0.002);
        simulator.SetSamplingTimestepMultiple(sample_rate); // Every 4 hours
        simulator.SetEndTime(end_time);

        // Pass an adaptive numerical method to the simulation
        boost::shared_ptr<AbstractNumericalMethod<2,2> > p_method(new ForwardEulerNumericalMethod<2,2>());
        p_method->SetUseAdaptiveTimestep(true);
        simulator.SetNumericalMethod(p_method);

        simulator.SetOutputDivisionLocations(true);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForceWithMinDistanceItem<2>, p_force);
        // MAKE_PTR(GeneralisedLinearSpringForceWithMinDistanceItem_mod<2>, p_force);
        // p_force->SetMeinekeSpringStiffness(30); //2.7 //15
        p_force->SetForceLawType(force_law);
        p_force->SetMeinekeSpringStiffness(spring_stiffness); //2.7 //15
        p_force->SetMeinekeDivisionRestingSpringLength(division_separation); //2.7
        p_force->SetMeinekeSpringGrowthDuration(2); //2.7
        p_force->SetCutOffLength(cut_off_length);
     
        simulator.AddForce(p_force);

        // MAKE_PTR(RandomMotionForce<2>, p_random_motion_force);
        // p_random_motion_force->SetMovementParameter(0.001);
        // simulator.AddForce(p_random_motion_force);


        // MAKE_PTR(ContactInhibitionModifier_mod<2>, p_contact_inhibition_modifier);
        MAKE_PTR(ContactInhibitionModifier<2>, p_contact_inhibition_modifier);
        simulator.AddSimulationModifier(p_contact_inhibition_modifier);

        // output the cell data
        // cell_population.AddCellWriter<CellDataItemWriter>();
        cell_population.AddPopulationWriter<NodeLocationWriter>();
        // cell_population.AddCellWriter<CellContactInhibitionWriter>();

        boost::shared_ptr<CellDataItemWriter<2,2> > p1_writer(new CellDataItemWriter<2,2>("FreeAreaFraction"));
        cell_population.AddCellWriter(p1_writer);

        boost::shared_ptr<CellDataItemWriter<2,2> > p2_writer(new CellDataItemWriter<2,2>("FreeSurfaceFraction"));
        cell_population.AddCellWriter(p2_writer);

        boost::shared_ptr<CellDataItemWriter<2,2> > p3_writer(new CellDataItemWriter<2,2>("growth inhibited"));
        cell_population.AddCellWriter(p3_writer);

        boost::shared_ptr<CellDataItemWriter<2,2> > p4_writer(new CellDataItemWriter<2,2>("Radius"));
        cell_population.AddCellWriter(p4_writer);

        boost::shared_ptr<CellDataItemWriter<2,2> > p5_writer(new CellDataItemWriter<2,2>("target area"));
        cell_population.AddCellWriter(p5_writer);

        boost::shared_ptr<CellDataItemWriter<2,2> > p6_writer(new CellDataItemWriter<2,2>("birth age"));
        cell_population.AddCellWriter(p6_writer);

        boost::shared_ptr<CellDataItemWriter<2,2> > p7_writer(new CellDataItemWriter<2,2>("Deformable Radius"));
        cell_population.AddCellWriter(p7_writer);

        simulator.Solve();

    }

    // void xTest2DMonolayerWithoutDiffusionSingleCell()
    // {

    //     double spring_stiffness = 0.0;
    //     // create a string identifier to take the string values of either linear, quadratic or log force law
    //     // std::string force_law = "quadratic"; // "log"; // "quadratic"; // "linear";
    //     TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-p_force_law"));
    //     std::string force_law = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-p_force_law");
    //     if (force_law == "linear")
    //     {
    //         spring_stiffness = 18.2816648;
    //     }
    //     else if (force_law == "quadratic")
    //     {
    //         // spring_stiffness = 88.54;
    //         spring_stiffness = 155.0;
    //     }
    //     else if (force_law == "log")
    //     {
    //         spring_stiffness = 15.6556;
    //     }
    //     PRINT_VARIABLE(force_law);
    //     PRINT_VARIABLE(spring_stiffness);

    //     // reset the random number generator
    //     double random_seed;
    //     try
    //     {
    //       TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-random_seed"));
    //       random_seed = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-random_seed");
    //     }
    //     catch(const std::exception& e)
    //     {
    //       random_seed = 0.0;
    //     }
    //     RandomNumberGenerator::Instance()->Reseed(random_seed);

    //     TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-p_beta"));
    //     double p_beta = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-p_beta");

    //     TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-p_gamma"));
    //     double p_gamma = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-p_gamma");

    //     TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-end_time"));
    //     double end_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-end_time");

    //     TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-sample_rate"));
    //     double sample_rate = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-sample_rate");

    //     TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-output_name"));
    //     std::string output_name = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-output_name");

    //     TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-cut_off_length"));
    //     double cut_off_length = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-cut_off_length");


    //     // static const double end_time = 200; //28*24; // 28 days first 14 days and second 14 days can be separated 


    //     NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        
    //     std::vector<double> cell_1 = {0.0, 0.0};
    //     Node<2> node1(0, cell_1.data(), false);
    //     p_mesh->AddNode(&node1);

    //     // std::vector<double> cell_2 = {0.5, 0.0};
    //     // Node<2> node2(1, cell_2.data(), false);
    //     // p_mesh->AddNode(&node2);

    //     // std::vector<double> cell_3 = {0.0, 0.5};
    //     // Node<2> node3(2, cell_3.data(), false);
    //     // p_mesh->AddNode(&node3);
      
    //     // std::vector<double> cell_4 = {-0.5, 0.0};
    //     // Node<2> node4(3, cell_4.data(), false);
    //     // p_mesh->AddNode(&node4);

    //     // std::vector<double> cell_5 = {0.0, -0.5};
    //     // Node<2> node5(4, cell_5.data(), false);
    //     // p_mesh->AddNode(&node5);

    //     p_mesh->SetMaximumInteractionDistance(1.1*cut_off_length);

    //     std::vector<CellPtr> cells;
    //     GenerateCells(p_mesh->GetNumNodes(),cells, false, p_beta, p_gamma);

    //     double division_separation = 2.0*(std::sqrt(1.0/M_PI));

    //     VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
    //     cell_population.AddCellWriter<CellIdWriter>();
    //     cell_population.AddCellWriter<CellAgesWriter>();
    //     cell_population.AddCellWriter<CellMutationStatesWriter>();
    //     cell_population.AddCellWriter<CellVolumesWriter>();
    //     // cell_population.AddCellWriter<CellContactInhibitionWriter>();
    //     cell_population.AddPopulationWriter<TissueWidthWriter>();
    //     cell_population.SetMeinekeDivisionSeparation(division_separation);
    //     cell_population.SetUseVariableRadii(true);

    //     OffLatticeSimulationWithPopulationBasedStoppingEvent simulator(cell_population);
    //     simulator.SetOutputDirectory(output_name);
    //     simulator.SetDt(0.002);
    //     simulator.SetSamplingTimestepMultiple(sample_rate); // Every 4 hours
    //     simulator.SetEndTime(end_time);

    //     // Pass an adaptive numerical method to the simulation
    //     boost::shared_ptr<AbstractNumericalMethod<2,2> > p_method(new ForwardEulerNumericalMethod<2,2>());
    //     p_method->SetUseAdaptiveTimestep(true);
    //     simulator.SetNumericalMethod(p_method);

    //     simulator.SetOutputDivisionLocations(true);

    //     // Create a force law and pass it to the simulation
    //     MAKE_PTR(GeneralisedLinearSpringForceWithMinDistanceItem<2>, p_force);
    //     // MAKE_PTR(GeneralisedLinearSpringForceWithMinDistanceItem_mod<2>, p_force);
    //     // p_force->SetMeinekeSpringStiffness(30); //2.7 //15
    //     p_force->SetForceLawType(force_law);
    //     p_force->SetMeinekeSpringStiffness(spring_stiffness); //2.7 //15
    //     p_force->SetMeinekeDivisionRestingSpringLength(division_separation); //2.7
    //     p_force->SetMeinekeSpringGrowthDuration(2); //2.7
    //     p_force->SetCutOffLength(cut_off_length);
     
    //     simulator.AddForce(p_force);

    //     // MAKE_PTR(ContactInhibitionModifier_mod<2>, p_contact_inhibition_modifier);
    //     MAKE_PTR(ContactInhibitionModifier<2>, p_contact_inhibition_modifier);
    //     simulator.AddSimulationModifier(p_contact_inhibition_modifier);

    //     // output the cell data
    //     // cell_population.AddCellWriter<CellDataItemWriter>();
    //     cell_population.AddPopulationWriter<NodeLocationWriter>();
    //     // cell_population.AddCellWriter<CellContactInhibitionWriter>();

    //     boost::shared_ptr<CellDataItemWriter<2,2> > p1_writer(new CellDataItemWriter<2,2>("FreeAreaFraction"));
    //     cell_population.AddCellWriter(p1_writer);

    //     boost::shared_ptr<CellDataItemWriter<2,2> > p2_writer(new CellDataItemWriter<2,2>("FreeSurfaceFraction"));
    //     cell_population.AddCellWriter(p2_writer);

    //     boost::shared_ptr<CellDataItemWriter<2,2> > p3_writer(new CellDataItemWriter<2,2>("growth inhibited"));
    //     cell_population.AddCellWriter(p3_writer);

    //     boost::shared_ptr<CellDataItemWriter<2,2> > p4_writer(new CellDataItemWriter<2,2>("Radius"));
    //     cell_population.AddCellWriter(p4_writer);

    //     boost::shared_ptr<CellDataItemWriter<2,2> > p5_writer(new CellDataItemWriter<2,2>("target area"));
    //     cell_population.AddCellWriter(p5_writer);

    //     boost::shared_ptr<CellDataItemWriter<2,2> > p6_writer(new CellDataItemWriter<2,2>("birth age"));
    //     cell_population.AddCellWriter(p6_writer);

    //     boost::shared_ptr<CellDataItemWriter<2,2> > p7_writer(new CellDataItemWriter<2,2>("Deformable Radius"));
    //     cell_population.AddCellWriter(p7_writer);

    //     simulator.Solve();

    // }

};

#endif /* TEST02MONOLAYERGROWTH_HPP_ */



