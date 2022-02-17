/*

Copyright (c) 2005-2019, University of Oxford.
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

#include "MeshRemeshModifier.hpp"
#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"
#include "VtkMeshWriter.hpp"

#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StromalCellMutationState.hpp"
#include "Debug.hpp"



template<unsigned DIM>
MeshRemeshModifier<DIM>::MeshRemeshModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mOutputDirectory(""),
    mWidth(0.0),
    mDepth(0.0)
{
}

template<unsigned DIM>
void MeshRemeshModifier<DIM>::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;
}

template<unsigned DIM>
void MeshRemeshModifier<DIM>::SetWidth(double width)
{
	mWidth = width;
}

template<unsigned DIM>
void MeshRemeshModifier<DIM>::SetDepth(double depth)
{
	mDepth = depth;
}

// std::string MeshModifier::GetOutputDirectory()
// {
// 	return mOutputDirectory;
// }

template<unsigned DIM>
MeshRemeshModifier<DIM>::~MeshRemeshModifier()
{
}

template<unsigned DIM>
void MeshRemeshModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
}

template<unsigned DIM>
void MeshRemeshModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
    std::vector<Node<DIM>*> extended_nodes_no_ghosts(3*rCellPopulation.GetNumRealCells());
    // std::vector<Node<DIM>*> extended_nodes_no_ghosts(rCellPopulation.GetNumRealCells());
    int count = 0;

    MutableMesh<DIM,DIM>* mpMeshNoGhosts;
    if (mpMeshNoGhosts != NULL)
    {
    	delete mpMeshNoGhosts;
    }

    std::map<unsigned, unsigned> MeshNodeIndexMap;

    // Dom - Create a copy of original mesh
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // First, create and store a copy of this real node and cell

        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
		// PRINT_2_VARIABLES(real_node_index,count);

        // Create a copy of the node corresponding to this cell and store it
        Node<DIM>* p_real_node = new Node<DIM>(real_node_index, real_node_location);
        extended_nodes_no_ghosts[count] = p_real_node;
        MeshNodeIndexMap[count] = real_node_index;
        count++;
        
    }

    // Dom - Create a copy of original mesh
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // First, create and store a copy of this real node and cell

        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        c_vector<double,DIM> image_node_location_width = real_node_location;
        if (real_node_location[0] >= mWidth*0.5)
        {
            image_node_location_width[0] -= mWidth;
        }
        else if (real_node_location[0] <  mWidth*0.5)
        {
            image_node_location_width[0] += mWidth;
        }

        Node<DIM>* p_real_node_width = new Node<DIM>(count, image_node_location_width);

        extended_nodes_no_ghosts[count] = p_real_node_width;
        MeshNodeIndexMap[count] = real_node_index;
        count++;


        c_vector<double,DIM> image_node_location_depth = real_node_location;
        if (real_node_location[1] >= mDepth*0.5)
        {
            image_node_location_depth[1] -= mDepth;
        }
        else if (real_node_location[1] <  mDepth*0.5)
        {
            image_node_location_depth[1] += mDepth;
        }

        Node<DIM>* p_real_node_depth = new Node<DIM>(count, image_node_location_depth);

        extended_nodes_no_ghosts[count] = p_real_node_depth;
        MeshNodeIndexMap[count] = real_node_index;
        count++;
        
        
    }

    mpMeshNoGhosts = new MutableMesh<DIM,DIM>(extended_nodes_no_ghosts);

    // PRINT_2_VARIABLES(rCellPopulation.GetNumRealCells() , mpMeshNoGhosts->GetNumNodes());
    
    MutableMesh<2,3>* pMutableMesh = new MutableMesh<2,3>();

    std::vector<Node<3>*> nodes;

    DomMeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);


    for (unsigned node_index=0; node_index<mpMeshNoGhosts->GetNumNodes(); node_index++)
    {
        const c_vector<double, 3>& r_location = mpMeshNoGhosts->GetNode(node_index)->rGetLocation();
        nodes.push_back(new Node<3>(node_index, r_location));
        pMutableMesh->AddNode(new Node<3>(node_index, r_location));
    }


    int element_iter = 0;
    double dist_threshold = 4.0;
    int num_nodes_original_mesh = rCellPopulation.GetNumRealCells();

    for (unsigned elem_index=0; elem_index<mpMeshNoGhosts->GetNumElements(); elem_index++) 
    { 
        bool element_contains_ghost_nodes = false;

        // Get a pointer to the element
        Element<DIM,DIM>* p_element = mpMeshNoGhosts->GetElement(elem_index);
        
        int Node0Index = p_element->GetNodeGlobalIndex(0);
        int Node1Index = p_element->GetNodeGlobalIndex(1);
        int Node2Index = p_element->GetNodeGlobalIndex(2);
        int Node3Index = p_element->GetNodeGlobalIndex(3);

        int node_intex[4] = {Node0Index, Node1Index, Node2Index, Node3Index};

        int node0GlobalIndex = MeshNodeIndexMap[Node0Index];
        int node1GlobalIndex = MeshNodeIndexMap[Node1Index];
        int node2GlobalIndex = MeshNodeIndexMap[Node2Index];
        int node3GlobalIndex = MeshNodeIndexMap[Node3Index];
        int node_global_intex[4] = {node0GlobalIndex, node1GlobalIndex, node2GlobalIndex, node3GlobalIndex};

        // PRINT_4_VARIABLES(Node0Index, Node1Index, Node2Index, Node3Index);
        // PRINT_2_VARIABLES(mpMeshNoGhosts->GetNumElements() , elem_index);

        // int number_of_ghosts = (p_tissue->IsGhostNode(node0GlobalIndex)) + (p_tissue->IsGhostNode(node1GlobalIndex)) + (p_tissue->IsGhostNode(node2GlobalIndex)) + (p_tissue->IsGhostNode(node3GlobalIndex));
        
        int part_of_original_mesh = 0;
        part_of_original_mesh = (Node0Index <= num_nodes_original_mesh) + (Node1Index <= num_nodes_original_mesh) + (Node2Index <= num_nodes_original_mesh) + (Node3Index <= num_nodes_original_mesh);
        // PRINT_VARIABLE(part_of_original_mesh);
        if (part_of_original_mesh > 1 )
        {
            for(int j=0; j<4; j++)
            {
                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_global_intex[(0+j)%4]);
                boost::shared_ptr<AbstractCellMutationState> p_state = p_cell->GetMutationState();

                bool is_0_stromal_cell = (p_state->IsType<StromalCellMutationState>()==true);

                if(is_0_stromal_cell) 
                {
                    CellPtr p_cell_1 = rCellPopulation.GetCellUsingLocationIndex(node_global_intex[(1+j)%4]);
                    boost::shared_ptr<AbstractCellMutationState> p_state_1 = p_cell_1->GetMutationState();

                    CellPtr p_cell_2 = rCellPopulation.GetCellUsingLocationIndex(node_global_intex[(2+j)%4]);
                    boost::shared_ptr<AbstractCellMutationState> p_state_2 = p_cell_2->GetMutationState();
                    
                    CellPtr p_cell_3 = rCellPopulation.GetCellUsingLocationIndex(node_global_intex[(3+j)%4]);
                    boost::shared_ptr<AbstractCellMutationState> p_state_3 = p_cell_3->GetMutationState();
                    
                    int number_of_epithelial_cell = (p_state_1->IsType<WildTypeCellMutationState>()==true) + (p_state_2->IsType<WildTypeCellMutationState>()==true) + (p_state_3->IsType<WildTypeCellMutationState>()==true);

                    if (number_of_epithelial_cell == 3)
                    {
                        // PRINT_3_VARIABLES(node_intex[(0+j)%4],node_intex[(1+j)%4],node_intex[(2+j)%4] );

                        const c_vector<double, 3>& r_location1 = rCellPopulation.GetNode(node_global_intex[(1+j)%4])->rGetLocation();
                        const c_vector<double, 3>& r_location2 = rCellPopulation.GetNode(node_global_intex[(2+j)%4])->rGetLocation();
                        const c_vector<double, 3>& r_location3 = rCellPopulation.GetNode(node_global_intex[(3+j)%4])->rGetLocation();

                        // c_vector<double, 3> difference12 = r_location1 - r_location2;
                        double distance_between_nodes_12 = norm_2(r_location1 - r_location2);

                        // c_vector<double, 3> difference13 = r_location1 - r_location3;
                        double distance_between_nodes_13 = norm_2(r_location1 - r_location3);

                        // c_vector<double, 3> difference23 = r_location2 - r_location3;
                        double distance_between_nodes_23 = norm_2(r_location2 - r_location3);

                        // We do this because the nodes may move after, due to periodicity, but the tesselation happens at the begining of the timestep
                        if(distance_between_nodes_12 <= dist_threshold && distance_between_nodes_13 <= dist_threshold && distance_between_nodes_23 <= dist_threshold)
                        {
                            std::vector<Node<3>*> ElementNodes0;
                            ElementNodes0.push_back(nodes[node_intex[(1+j)%4]]);
                            ElementNodes0.push_back(nodes[node_intex[(2+j)%4]]);
                            ElementNodes0.push_back(nodes[node_intex[(3+j)%4]]);
                            pMutableMesh->AddElement(new Element<2,3>(element_iter,ElementNodes0));
                            element_iter++;

                        }
                    }
                    
                }
                
                
            }
        }

    }


    std::ostringstream time_string;
    time_string << SimulationTime::Instance()->GetTimeStepsElapsed();
    std::string results_file = "mono_mesh_results_" + time_string.str();
    VtkMeshWriter<2,3>* p_vtk_mesh_writer = new VtkMeshWriter<2,3>(mOutputDirectory, results_file, false);
    p_vtk_mesh_writer->WriteFilesUsingMesh(*pMutableMesh);

    delete pMutableMesh;
    delete mpMeshNoGhosts;

    // delete pOriginalMesh;
}



template<unsigned DIM>
void MeshRemeshModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateAtEndOfOutputTimeStep(rCellPopulation);


}


template<unsigned DIM>
void MeshRemeshModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
// template class MeshModifier<1>;
// template class MeshModifier<2>;
template class MeshRemeshModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeshRemeshModifier)

