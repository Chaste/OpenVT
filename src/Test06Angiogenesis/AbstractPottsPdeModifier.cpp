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

#include "AbstractPottsPdeModifier.hpp"
#include "ReplicatableVector.hpp"
#include "LinearBasisFunction.hpp"
#include "Debug.hpp"
#include "VtkMeshWriter.hpp"
#include "PottsMesh.hpp"

template<unsigned DIM>
AbstractPottsPdeModifier<DIM>::AbstractPottsPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde,
                                                                boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition,
                                                                bool isNeumannBoundaryCondition,
                                                                boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid,
                                                                double stepSize,
                                                                Vec solution)
    : AbstractPdeModifier<DIM>(pPde,
                               pBoundaryCondition,
                               isNeumannBoundaryCondition,
                               solution),
      mpMeshCuboid(pMeshCuboid),
      mStepSize(stepSize),
      mSetBcsOnBoxBoundary(true)
{
    if (pMeshCuboid)
    {
        // We only need to generate mpFeMesh once, as it does not vary with time
        this->GenerateFeMesh(mpMeshCuboid, mStepSize);
        this->mDeleteFeMesh = true;
        //initialise the boundary nodes
        this->mIsDirichletBoundaryNode = std::vector<double>(this->mpFeMesh->GetNumNodes(), 0.0);
    }
}

template<unsigned DIM>
AbstractPottsPdeModifier<DIM>::~AbstractPottsPdeModifier()
{
}

template<unsigned DIM>
double AbstractPottsPdeModifier<DIM>::GetStepSize()
{
     return mStepSize;
}

template<unsigned DIM>
void AbstractPottsPdeModifier<DIM>::SetBcsOnBoxBoundary(bool setBcsOnBoxBoundary)
{
    mSetBcsOnBoxBoundary = setBcsOnBoxBoundary;
}

template<unsigned DIM>
bool AbstractPottsPdeModifier<DIM>::AreBcsSetOnBoxBoundary()
{
    return mSetBcsOnBoxBoundary;
}

template<unsigned DIM>
void AbstractPottsPdeModifier<DIM>::ConstructBoundaryConditionsContainerHelper(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                                                                                   std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > pBcc)
{
    // Check both meshes are the same
        if (this->mpFeMesh->GetNumNodes() != rCellPopulation.GetNumNodes())
        {
            PRINT_2_VARIABLES(this->mpFeMesh->GetNumNodes(), rCellPopulation.GetNumNodes());

            NEVER_REACHED;
        }

    // Reset mIsDirichletBoundaryNode
    for (unsigned i=0; i<this->mpFeMesh->GetNumNodes(); i++)
    {
       this->mIsDirichletBoundaryNode[i] = 0.0;
    }
// TODO make this use PottsElements to define the tissue.
    if (!this->mSetBcsOnBoxBoundary)
    {
        NEVER_REACHED;
        // Set pde nodes as boundary node if elements dont contain cells or nodes aren't within a specified distance of a cell centre
        
        // // Get the set of coarse element indices that contain cells
        // std::set<unsigned> coarse_element_indices_in_map;
        // for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
        //     cell_iter != rCellPopulation.End();
        //     ++cell_iter)
        // {
        //     coarse_element_indices_in_map.insert(this->mCellPdeElementMap[*cell_iter]);
        // }

        // // Find the node indices associated with elements whose indices are NOT in the set coarse_element_indices_in_map
        // std::set<unsigned> coarse_mesh_boundary_node_indices;
        // for (unsigned i=0; i<this->mpFeMesh->GetNumElements(); i++)
        // {
        //     if (coarse_element_indices_in_map.find(i) == coarse_element_indices_in_map.end())
        //     {
        //         Element<DIM,DIM>* p_element = this->mpFeMesh->GetElement(i);
        //         for (unsigned j=0; j<DIM+1; j++)
        //         {
        //             unsigned node_index = p_element->GetNodeGlobalIndex(j);
        //             coarse_mesh_boundary_node_indices.insert(node_index);
        //         }
        //     }
        // }

        // // Also remove nodes that are within the typical cell radius from the centre of a cell.
        // std::set<unsigned> nearby_node_indices;
        // for (std::set<unsigned>::iterator node_iter = coarse_mesh_boundary_node_indices.begin();
        //     node_iter != coarse_mesh_boundary_node_indices.end();
        //     ++node_iter)
        // {
        //     bool remove_node = false;

        //     c_vector<double,DIM> node_location = this->mpFeMesh->GetNode(*node_iter)->rGetLocation();

        //     for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
        //     cell_iter != rCellPopulation.End();
        //     ++cell_iter)
        //     {
        //         const ChastePoint<DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        //         double separation = norm_2(node_location - r_position_of_cell.rGetLocation());

        //         if (separation <= mTypicalCellRadius)
        //         {
        //             remove_node = true;
        //             break;
        //         }                
        //     }

        //     if (remove_node)
        //     {
        //         // Node near cell so set it to be removed from boundary set
        //         nearby_node_indices.insert(*node_iter);
        //     }
        // }

        // // Remove nodes that are near cells from boundary set
        // for (std::set<unsigned>::iterator node_iter = nearby_node_indices.begin();
        //     node_iter != nearby_node_indices.end();
        //     ++node_iter)
        // {
        //     coarse_mesh_boundary_node_indices.erase(*node_iter);
        // }


        // // Apply boundary condition to the nodes in the set coarse_mesh_boundary_node_indices
        // if (this->IsNeumannBoundaryCondition())
        // {
        //     // Neumann BSC not implemented yet for this geometry
        //     NEVER_REACHED;
        // }
        // else
        // {
        //     // Impose any Dirichlet boundary conditions
        //     for (std::set<unsigned>::iterator iter = coarse_mesh_boundary_node_indices.begin();
        //         iter != coarse_mesh_boundary_node_indices.end();
        //         ++iter)
        //     {
        //         pBcc->AddDirichletBoundaryCondition(this->mpFeMesh->GetNode(*iter), this->mpBoundaryCondition.get(), 0, false);
        //         this->mIsDirichletBoundaryNode[*iter] = 1.0;
        //     }
        // }
        
    }
    else // Apply BC at boundary of box domain FE mesh
    {
        if (this->IsNeumannBoundaryCondition())
        {
            // Impose any Neumann boundary conditions
            for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();
                 elem_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd();
                 ++elem_iter)
            {
                pBcc->AddNeumannBoundaryCondition(*elem_iter, this->mpBoundaryCondition.get());
            }
        }
        else
        {
            // Impose any Dirichlet boundary conditions
            for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
                 node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
                 ++node_iter)
            {
                pBcc->AddDirichletBoundaryCondition(*node_iter, this->mpBoundaryCondition.get());
                this->mIsDirichletBoundaryNode[(*node_iter)->GetIndex()] = 1.0;
            }
        }
    }
}



template<unsigned DIM>
void AbstractPottsPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractPdeModifier<DIM>::SetupSolve(rCellPopulation, outputDirectory);

    // InitialiseCellPdeElementMap(rCellPopulation);
}

template<unsigned DIM>
void AbstractPottsPdeModifier<DIM>::GenerateFeMesh(boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid, double stepSize)
{
    // Create a regular coarse tetrahedral mesh
    this->mpFeMesh = new TetrahedralMesh<DIM,DIM>();
//TODO change to use the PottsMesh
    GenerateAndReturnFeMesh(pMeshCuboid, stepSize, this->mpFeMesh);
}

template<unsigned DIM>
void AbstractPottsPdeModifier<DIM>::GenerateAndReturnFeMesh(boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid, double stepSize, TetrahedralMesh<DIM,DIM>* pMesh)
{
//TODO change to use the PottsMesh

    // Create a regular coarse tetrahedral mesh
    switch (DIM)
    {
        case 1:
            pMesh->ConstructRegularSlabMesh(stepSize, pMeshCuboid->GetWidth(0));
            break;
        case 2:
            pMesh->ConstructRegularSlabMesh(stepSize, pMeshCuboid->GetWidth(0), pMeshCuboid->GetWidth(1));
            break;
        case 3:
            pMesh->ConstructRegularSlabMesh(stepSize, pMeshCuboid->GetWidth(0), pMeshCuboid->GetWidth(1), pMeshCuboid->GetWidth(2));
            break;
        default:
            NEVER_REACHED;
    }

    // Get centroid of meshCuboid
    ChastePoint<DIM> upper = pMeshCuboid->rGetUpperCorner();
    ChastePoint<DIM> lower = pMeshCuboid->rGetLowerCorner();
    c_vector<double,DIM> centre_of_cuboid = 0.5*(upper.rGetLocation() + lower.rGetLocation());

    // Find the centre of the PDE mesh
    c_vector<double,DIM> centre_of_coarse_mesh = zero_vector<double>(DIM);
    for (unsigned i=0; i<pMesh->GetNumNodes(); i++)
    {
        centre_of_coarse_mesh += pMesh->GetNode(i)->rGetLocation();
    }
    centre_of_coarse_mesh /= pMesh->GetNumNodes();

    // Now move the mesh to the correct location
    pMesh->Translate(centre_of_cuboid - centre_of_coarse_mesh);
}

template<unsigned DIM>
void AbstractPottsPdeModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    PottsMesh<DIM>* p_static_cast_mesh = static_cast<PottsMesh<DIM>*>(&(rCellPopulation.rGetMesh()));
    
    // Store the PDE solution in an accessible form
    ReplicatableVector solution_repl(this->mSolution);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {
        // The cells are not nodes of the mesh, so we must loop over all contained nodes 
        double solution_at_cell = 0.0;
    
        unsigned cell_location_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        
        PottsElement<DIM>* p_element = p_static_cast_mesh->GetElement(cell_location_index);
        
        unsigned num_nodes_in_element = p_element->GetNumNodes();

        for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
        {
            // Find location of current node and add it to the centroid
            solution_at_cell += solution_repl[p_element->GetNodeGlobalIndex(local_index)];
        }

        solution_at_cell = solution_at_cell/num_nodes_in_element;

        cell_iter->GetCellData()->SetItem(this->mDependentVariableName, solution_at_cell);
    }   

    if (this->mOutputGradient)
    {
        // See other PDE Modifiers for how to implement. 
        NEVER_REACHED;
    }
}

template<unsigned DIM>
void AbstractPottsPdeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractPdeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class AbstractPottsPdeModifier<1>;
template class AbstractPottsPdeModifier<2>;
template class AbstractPottsPdeModifier<3>;
